#!/usr/bin/env bash
# Launch the staging EC2 (eu-west-1, default VPC): key pair, security group with SSH
# from the caller's current public IP only, and an m7g.xlarge with a 300 GB gp3 root.
# Everything is tagged Project=cnv_cfs so teardown_staging_ec2.sh can find it.
#
# Usage: bash scripts/01_staging/launch_staging_ec2.sh [--profile <AWS_PROFILE>]
set -euo pipefail

PROFILE=<AWS_PROFILE>
REGION=eu-west-1
NAME=cnv-cfs-staging
TYPE=m7g.xlarge
DISK_GB=300
KEY_FILE="$HOME/.ssh/${NAME}.pem"
EXPECT_ACCOUNT=<AWS_ACCOUNT_ID>
[[ "${1:-}" == "--profile" ]] && PROFILE="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
AWS=(aws --profile "$PROFILE" --region "$REGION")
TAGS="ResourceType=%s,Tags=[{Key=Name,Value=$NAME},{Key=Project,Value=cnv_cfs}]"

ACCOUNT=$("${AWS[@]}" sts get-caller-identity --query Account --output text)
[[ "$ACCOUNT" == "$EXPECT_ACCOUNT" ]] || { echo "account mismatch: $ACCOUNT" >&2; exit 1; }

# key pair (private key kept locally only)
if [[ ! -f "$KEY_FILE" ]]; then
  "${AWS[@]}" ec2 create-key-pair --key-name "$NAME" --key-type ed25519 \
    --tag-specifications "$(printf "$TAGS" key-pair)" \
    --query KeyMaterial --output text > "$KEY_FILE"
  chmod 600 "$KEY_FILE"
fi

# security group: SSH from this IP only
VPC=$("${AWS[@]}" ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
SG=$("${AWS[@]}" ec2 describe-security-groups --filters Name=group-name,Values="$NAME" \
     --query 'SecurityGroups[0].GroupId' --output text)
if [[ "$SG" == "None" ]]; then
  SG=$("${AWS[@]}" ec2 create-security-group --group-name "$NAME" --vpc-id "$VPC" \
       --description "SSH for cnv_cfs SRA staging" \
       --tag-specifications "$(printf "$TAGS" security-group)" --query GroupId --output text)
fi
MYIP=$(curl -sSf https://checkip.amazonaws.com | tr -d '[:space:]')
"${AWS[@]}" ec2 authorize-security-group-ingress --group-id "$SG" \
  --protocol tcp --port 22 --cidr "$MYIP/32" 2>/dev/null || true

AMI=$("${AWS[@]}" ssm get-parameter \
      --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
      --query Parameter.Value --output text)

IID=$("${AWS[@]}" ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" \
      --key-name "$NAME" --security-group-ids "$SG" \
      --block-device-mappings "DeviceName=/dev/xvda,Ebs={VolumeSize=$DISK_GB,VolumeType=gp3,DeleteOnTermination=true}" \
      --instance-initiated-shutdown-behavior terminate \
      --user-data "file://$HERE/ec2_userdata.sh" \
      --tag-specifications "$(printf "$TAGS" instance)" "$(printf "$TAGS" volume)" \
      --query 'Instances[0].InstanceId' --output text)
echo "instance: $IID  (waiting for running state)"
"${AWS[@]}" ec2 wait instance-running --instance-ids "$IID"
IP=$("${AWS[@]}" ec2 describe-instances --instance-ids "$IID" \
     --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
echo "public IP: $IP"
echo "ssh -i $KEY_FILE ec2-user@$IP"
