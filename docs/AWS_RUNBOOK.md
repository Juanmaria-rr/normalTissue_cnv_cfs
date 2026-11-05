# AWS runbook — staging PRJNA1003661 SRA files for HealthOmics

**Technical report and operating instructions.** First written 2026-10-23, from the pilot
staging run. This document lets the procedure be repeated by hand, without the helper
scripts and without any AI assistant. Each step gives the command, the expected result and
what to do if it fails.

The scripts in `scripts/01_staging/` automate exactly these steps. Where a script exists,
it is named, but the raw commands are always given.

---

## 0. What this procedure does, and why it is needed

| Fact | Consequence |
|---|---|
| The SRA run files live in the public bucket `s3://sra-pub-run-odp/sra/<SRR>/<SRR>`, region **us-east-1** | They are free to read, but in the wrong region for us |
| Our HealthOmics account runs in **eu-west-1**; an organisation SCP denies HealthOmics in us-east-1 | Workflows can only read inputs from an eu-west-1 bucket |
| Signed reads of `sra-pub-run-odp` with the work account return `AccessDenied` (org policy). **Unsigned** reads (`--no-sign-request`) work | Downloads must be unsigned; uploads must be signed. No server-side cross-region copy is possible |
| The CNIO workstation was unavailable, and the laptop uplink is ~12.5 MB/s | We use a temporary EC2 in eu-west-1 as the transfer host |
| The EC2 has no IAM role (we cannot create or list instance profiles), and exported laptop credentials expire every ~10–15 min | The instance authenticates with `aws login`, which gives a ~70 h self-refreshing session |

Data path: `sra-pub-run-odp` (us-east-1) → EC2 local disk (eu-west-1) →
`s3://<PROJECT_BUCKET>/cnv_cfs/sra/` (eu-west-1).
Egress from the open-data bucket is not charged, and ingress to S3 is free. The cost is the
instance, its disk and S3 storage (§9).

---

## 1. Accounts and identities

| | value |
|---|---|
| AWS account | **<AWS_ACCOUNT_ID>** |
| IAM user | `<IAM_USER>` |
| CLI profile | `<AWS_PROFILE>` (the profile `<IAM_USER>` reaches the same user) |
| Region | `eu-west-1` |
| Bucket / prefix | `s3://<PROJECT_BUCKET>/cnv_cfs/` |

**Do not use** the laptop's `default` profile (account <OTHER_AWS_ACCOUNT_ID>, no S3/omics rights). Do
not use the `tailor` profile either: it belongs to another organisation.

```bash
aws --version                         # must be aws-cli/2.30 or newer (aws login)
aws login --profile <AWS_PROFILE>  # browser sign-in + MFA
export AWS_PROFILE=<AWS_PROFILE> AWS_REGION=eu-west-1
aws sts get-caller-identity
```

Expected: `"Account": "<AWS_ACCOUNT_ID>"`, ARN ending in `user/<IAM_USER>`.

---

## 2. Inputs: which runs to copy

The list of runs comes from the metadata step (see `README.md`, *Reproduce*):

```bash
python3 scripts/00_metadata/build_cell_metadata.py --out-dir data/metadata
python3 scripts/00_metadata/select_pilot_cells.py --cells data/metadata/cells.tsv \
    --donors P08 P07 P36 P20 P18 --out-dir data/pilot
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{print $h["run_accession"]}' \
    data/pilot/pilot_runs.tsv | sort -u > data/pilot/runs_pilot.txt
wc -l data/pilot/runs_pilot.txt       # pilot: 12542
head -3 data/pilot/runs_pilot.txt > data/pilot/runs_test3.txt
```

One line per SRR accession. A cell can have two runs (top-up sequencing). Both are staged;
they are merged per cell later, in the pipeline.

Sanity check from the laptop that the source is readable unsigned:

```bash
aws s3 ls s3://sra-pub-run-odp/sra/SRR30459505/ --no-sign-request
# 2024-10-19 12:03:34    4917796 SRR30459505
```

---

## 3. Pre-flight permission checks (optional, read-only)

Run these before launching anything, to confirm that EC2 is allowed. A `--dry-run` that
answers `DryRunOperation: Request would have succeeded` means "allowed"; nothing is created.

```bash
AMI=$(aws ssm get-parameter \
      --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
      --query Parameter.Value --output text)
aws ec2 run-instances --dry-run --image-id "$AMI" --instance-type t4g.small
aws ec2 create-key-pair --dry-run --key-name dryrun-test
VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
aws ec2 create-security-group --dry-run --group-name dryrun-test --description x --vpc-id "$VPC"
```

Verified 2026-10-23: all three succeed. `iam:ListInstanceProfiles` and
`iam:SimulatePrincipalPolicy` are denied, so there is no usable instance role (§5).

---

## 4. Launch the transfer instance

Script: `bash scripts/01_staging/launch_staging_ec2.sh`. By hand:

```bash
NAME=cnv-cfs-staging
TAGS="ResourceType=%s,Tags=[{Key=Name,Value=$NAME},{Key=Project,Value=cnv_cfs}]"

# 4.1 key pair (private key stays on the laptop)
aws ec2 create-key-pair --key-name $NAME --key-type ed25519 \
    --tag-specifications "$(printf "$TAGS" key-pair)" \
    --query KeyMaterial --output text > ~/.ssh/$NAME.pem
chmod 600 ~/.ssh/$NAME.pem

# 4.2 security group: SSH from this machine's public IP only
VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
SG=$(aws ec2 create-security-group --group-name $NAME --vpc-id "$VPC" \
     --description "SSH for cnv_cfs SRA staging" \
     --tag-specifications "$(printf "$TAGS" security-group)" --query GroupId --output text)
MYIP=$(curl -s https://checkip.amazonaws.com)
aws ec2 authorize-security-group-ingress --group-id "$SG" --protocol tcp --port 22 --cidr "$MYIP/32"

# 4.3 instance: Graviton m7g.xlarge (4 vCPU, up to 12.5 Gbit/s), 300 GB gp3 root
AMI=$(aws ssm get-parameter \
      --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
      --query Parameter.Value --output text)
IID=$(aws ec2 run-instances --image-id "$AMI" --instance-type m7g.xlarge \
      --key-name $NAME --security-group-ids "$SG" \
      --block-device-mappings 'DeviceName=/dev/xvda,Ebs={VolumeSize=300,VolumeType=gp3,DeleteOnTermination=true}' \
      --instance-initiated-shutdown-behavior terminate \
      --user-data file://scripts/01_staging/ec2_userdata.sh \
      --tag-specifications "$(printf "$TAGS" instance)" "$(printf "$TAGS" volume)" \
      --query 'Instances[0].InstanceId' --output text)
aws ec2 wait instance-running --instance-ids "$IID"
IP=$(aws ec2 describe-instances --instance-ids "$IID" \
     --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
echo "ssh -i ~/.ssh/$NAME.pem ec2-user@$IP"
```

`ec2_userdata.sh` installs `tmux` and the **latest AWS CLI v2**. The Amazon Linux package is
too old for `aws login`. Wait until setup has finished (1–2 min):

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP 'test -f /var/tmp/userdata_done && aws --version'
```

Expected: `aws-cli/2.3x ...`. If the file does not exist yet, wait and retry.

**Sizing.** 300 GB holds the whole pilot (143 GB of `.sra` measured). The disk must hold
everything staged in one batch, because files are downloaded first and synced afterwards.

**Disk throughput.** A gp3 volume defaults to **125 MB/s**. That was the bottleneck of the
pilot's sha256 and upload phases: each reads the whole batch from disk, ~18 min per pass for
143 GB. For larger batches, add `Throughput=500` (MB/s, up to 1000) and e.g. `Iops=6000` to
the `Ebs={...}` mapping. The extra cost is small and only lasts while the volume exists.

---

## 5. Authenticate on the instance (the non-obvious step)

The instance has no role, so it signs in as your IAM user with `aws login`. Three details,
each of which failed once:

**5.1 Preset the region.** Otherwise `aws login` stops at an interactive "AWS Region"
prompt:

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP \
  'aws configure set region eu-west-1 --profile <AWS_PROFILE>'
```

**5.2 Run the login from a real terminal** (Terminal.app, iTerm), not from a tool that has no
TTY. It needs keyboard input:

```bash
ssh -t -i ~/.ssh/$NAME.pem ec2-user@$IP 'aws login --remote --profile <AWS_PROFILE>'
```

Open the printed URL in your browser, then sign in and complete MFA.

**5.3 Deliver the redirect to the instance.** After MFA the browser is redirected to
`http://127.0.0.1:<port>/...?code=...`, which fails with `ERR_CONNECTION_REFUSED`. That address
is your laptop, but the CLI waiting for it is on the instance. **Leave the login command
running**, copy the full URL from the address bar, and in a second terminal run:

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP "curl -s 'http://127.0.0.1:<port>/...<full URL>...'"
```

The reply is an HTML page containing *"Your credentials have been shared successfully"*, and
the first terminal finishes. Do it within a few minutes; the code expires.

(Alternative: if you can read the port from the `redirect_uri=` parameter of the login URL
*before* opening it, open a tunnel first with
`ssh -i ~/.ssh/$NAME.pem -N -L <port>:127.0.0.1:<port> ec2-user@$IP`. The redirect then
reaches the instance by itself.)

**5.4 Attach the session to the named profile.** In our run the session was written to the
`default` profile, even though `--profile` was given. Attach it:

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP '
  aws configure set login_session <IAM_USER_ARN> --profile <AWS_PROFILE>
  aws sts get-caller-identity --profile <AWS_PROFILE> --query Arn --output text
  aws s3 ls s3://<PROJECT_BUCKET>/ --profile <AWS_PROFILE>'
```

Expected: the user ARN, then the bucket's top-level prefixes. The session refreshes itself for
about 70 h. After that, repeat §5.2–5.3.

---

## 6. Copy the run files

Script: `scripts/01_staging/stage_sra_to_s3.sh` (copy it and the run list to the instance).

```bash
scp -i ~/.ssh/$NAME.pem scripts/01_staging/stage_sra_to_s3.sh \
    data/pilot/runs_pilot.txt data/pilot/runs_test3.txt ec2-user@$IP:/data/
```

### 6.1 What the script does, as raw commands

For each SRR, in parallel (`xargs -P 64`):

```bash
# expected size at the source (unsigned)
aws s3api head-object --no-sign-request --region us-east-1 \
    --bucket sra-pub-run-odp --key sra/$RUN/$RUN --query ContentLength --output text
# download (unsigned); skip if a local file of that size already exists
aws s3 cp --no-sign-request --region us-east-1 \
    s3://sra-pub-run-odp/sra/$RUN/$RUN /data/sra/$RUN/$RUN
# compare: stat -c %s /data/sra/$RUN/$RUN  ==  expected size
```

Each run is recorded as `OK`, `SKIP` (already present), `SIZE_MISMATCH`, `DOWNLOAD_FAILED` or
`MISSING_AT_SOURCE` in `/data/logs/stage_<timestamp>.download.tsv`. Only `OK`/`SKIP` runs go on.

Then:

```bash
# manifest: run, bytes, sha256
sha256sum /data/sra/$RUN/$RUN          # per run, collected into stage_<ts>.manifest.tsv
# upload (signed)
aws s3 sync /data/sra s3://<PROJECT_BUCKET>/cnv_cfs/sra --profile <AWS_PROFILE>
aws s3 cp /data/logs/stage_<ts>.manifest.tsv \
    s3://<PROJECT_BUCKET>/cnv_cfs/sra/_manifest/ --profile <AWS_PROFILE>
# verify: every (run, size) pair of the manifest exists at the destination
aws s3 ls s3://<PROJECT_BUCKET>/cnv_cfs/sra/ --recursive --profile <AWS_PROFILE>
```

The last line printed is `STAGING VERIFIED` (all manifest runs present with the same size) or
`STAGING INCOMPLETE` (exit code 1).

### 6.2 Test with three runs first

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP 'cd /data && bash stage_sra_to_s3.sh \
  --runs runs_test3.txt --dest s3://<PROJECT_BUCKET>/cnv_cfs/_test/sra \
  --profile <AWS_PROFILE> --workdir /data/test/sra --jobs 4'
# expected: "3 OK" ... "missing/size-mismatch: 0" ... "STAGING VERIFIED"
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP '
  aws s3 rm --recursive s3://<PROJECT_BUCKET>/cnv_cfs/_test/ --profile <AWS_PROFILE>
  rm -rf /data/test'
```

### 6.3 The real run, inside tmux

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP
cd /data
tmux new -s staging
bash stage_sra_to_s3.sh --runs runs_pilot.txt \
     --dest s3://<PROJECT_BUCKET>/cnv_cfs/sra \
     --profile <AWS_PROFILE> --workdir /data/sra --jobs 64 2>&1 | tee /data/staging_pilot.log
# Ctrl-b d to detach; `tmux attach -t staging` to return
```

The transfer survives closing the laptop or dropping SSH.

### 6.4 Watching progress

```bash
ssh -i ~/.ssh/$NAME.pem ec2-user@$IP '
  find /data/sra -type f | wc -l          # runs downloaded so far
  du -sh /data/sra
  tail -3 /data/staging_pilot.log
  tmux has-session -t staging && echo running || echo finished'
```

Measured in the pilot: **~2.3–2.5 runs/s, ~46 MB/s**, about 1.5 h for 12,542 runs. The
bottleneck is CPU (two AWS CLI processes start per run; load average ~42 on 4 vCPU), not
network. For the full dataset (~99k runs, ~1.4 TB), use a larger instance (e.g. `m7g.4xlarge`
with `--jobs 128`) and a bigger disk, or stage in batches of a few donors.

### 6.5 Re-running after an interruption

The script is idempotent. Run the same command again: runs already on disk with the correct
size are `SKIP`ped, and `aws s3 sync` uploads only what is missing. If credentials have expired
(`NoCredentials` / `ExpiredToken` in the log), redo §5.2–5.4 first.

---

## 7. Final verification (independent of the script)

From the laptop:

```bash
aws s3 ls s3://<PROJECT_BUCKET>/cnv_cfs/sra/ --recursive --profile <AWS_PROFILE> \
  | awk '{n=split($4,p,"/"); if (p[n]==p[n-1]) print p[n]}' | sort > /tmp/at_dest.txt
comm -23 <(sort -u data/pilot/runs_pilot.txt) /tmp/at_dest.txt > /tmp/not_staged.txt
wc -l /tmp/at_dest.txt /tmp/not_staged.txt
```

Every entry of `not_staged.txt` must appear as `MISSING_AT_SOURCE`, `DOWNLOAD_FAILED` or
`SIZE_MISMATCH` in the instance's `download.tsv`. A run that appears in neither file means
something went wrong; investigate before tearing down. Copy the logs back for the record:

```bash
scp -i ~/.ssh/$NAME.pem -r ec2-user@$IP:/data/logs data/pilot/staging_logs
```

---

## 8. Tear down (do not skip: the instance bills by the hour)

Script: `bash scripts/01_staging/teardown_staging_ec2.sh`. By hand:

```bash
IIDS=$(aws ec2 describe-instances \
       --filters Name=tag:Name,Values=$NAME Name=instance-state-name,Values=pending,running,stopping,stopped \
       --query 'Reservations[].Instances[].InstanceId' --output text)
aws ec2 terminate-instances --instance-ids $IIDS
aws ec2 wait instance-terminated --instance-ids $IIDS
SG=$(aws ec2 describe-security-groups --filters Name=group-name,Values=$NAME --query 'SecurityGroups[0].GroupId' --output text)
aws ec2 delete-security-group --group-id "$SG"
aws ec2 delete-key-pair --key-name $NAME
rm -f ~/.ssh/$NAME.pem
```

The root volume is deleted with the instance (`DeleteOnTermination=true`). Data already in S3 is
not affected. Check that nothing is left:

```bash
aws ec2 describe-instances --filters Name=tag:Project,Values=cnv_cfs \
    --query 'Reservations[].Instances[].[InstanceId,State.Name]' --output text   # only "terminated"
aws ec2 describe-volumes --filters Name=tag:Project,Values=cnv_cfs --output text  # empty
```

---

## 9. Cost (eu-west-1 list prices, approximate)

| item | rate | pilot |
|---|---|---|
| m7g.xlarge on-demand | ~0.18 $/h | ~2 h → ~0.4 $ |
| gp3 300 GB | ~0.088 $/GB-month | a few hours → < 0.1 $ |
| Transfer out of the open-data bucket | not charged | 0 $ |
| Transfer into S3 | free | 0 $ |
| S3 Standard storage | ~0.023 $/GB-month | ~215 GB → ~5 $/month while kept |

Delete or tier the `.sra` copies when the project no longer needs them.

---

## 10. Failure modes

| symptom | meaning | action |
|---|---|---|
| `AccessDenied` listing `sra-pub-run-odp` | you sent a **signed** request; the org policy blocks it | add `--no-sign-request` |
| `aws login` asks "AWS Region [us-east-1]" | no region configured on the instance | §5.1 |
| `Input is not a terminal` / `Pseudo-terminal will not be allocated` | login launched without a TTY | run it from a real terminal, with `ssh -t` |
| browser: `127.0.0.1 ERR_CONNECTION_REFUSED` after MFA | redirect went to the laptop, not to the instance | §5.3 (`curl` the URL on the instance) |
| `NoCredentials` right after a successful login | session saved under `default` | §5.4 |
| `NoCredentials` / `ExpiredToken` mid-transfer | ~70 h session ended | §5.2–5.4, then re-run (§6.5) |
| `MISSING_AT_SOURCE` for some runs | run not present in the open-data bucket (cold or withdrawn) | list them, then fetch with `prefetch`/`fasterq-dump` from NCBI or drop them; record the decision |
| `ssh: connection timed out` | your public IP changed (VPN, network) | `aws ec2 authorize-security-group-ingress ... --cidr <new IP>/32` |
| disk full | batch larger than the disk | bigger volume, or stage fewer donors per batch |
| `iam:ListInstanceProfiles` / `SimulatePrincipalPolicy` denied | the user has no IAM read rights | expected; this is why §5 exists |

---

## 11. What this procedure must never do

- Never write to `s3://<OTHER_BUCKET>/cnio_cegen/`. It is the Deep Archive backup of
  another project.
- Never put credentials in user data, in the repository or in S3. The only credential on the
  instance is the `aws login` session, and it is destroyed with the instance.
- Never leave the instance running after verification (§8).
