#!/usr/bin/env bash
# Terminate the staging instance and delete its security group and key pair.
# Only touches resources tagged/named cnv-cfs-staging. Data already in S3 is not touched.
#
# Usage: bash scripts/01_staging/teardown_staging_ec2.sh [--profile <AWS_PROFILE>]
set -euo pipefail

PROFILE=<AWS_PROFILE>
REGION=eu-west-1
NAME=cnv-cfs-staging
[[ "${1:-}" == "--profile" ]] && PROFILE="$2"
AWS=(aws --profile "$PROFILE" --region "$REGION")

IIDS=$("${AWS[@]}" ec2 describe-instances \
       --filters Name=tag:Name,Values="$NAME" Name=tag:Project,Values=cnv_cfs \
                 Name=instance-state-name,Values=pending,running,stopping,stopped \
       --query 'Reservations[].Instances[].InstanceId' --output text)
if [[ -n "$IIDS" ]]; then
  echo "terminating: $IIDS"
  "${AWS[@]}" ec2 terminate-instances --instance-ids $IIDS > /dev/null
  "${AWS[@]}" ec2 wait instance-terminated --instance-ids $IIDS
fi
SG=$("${AWS[@]}" ec2 describe-security-groups --filters Name=group-name,Values="$NAME" \
     --query 'SecurityGroups[0].GroupId' --output text)
[[ "$SG" != "None" ]] && "${AWS[@]}" ec2 delete-security-group --group-id "$SG" && echo "deleted SG $SG"
"${AWS[@]}" ec2 delete-key-pair --key-name "$NAME" && echo "deleted key pair $NAME"
rm -f "$HOME/.ssh/${NAME}.pem"
echo "teardown complete"
