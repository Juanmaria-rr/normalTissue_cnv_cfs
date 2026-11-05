#!/bin/bash
# User data for the staging instance (Amazon Linux 2023, arm64): latest AWS CLI v2
# (`aws login` needs >= 2.30; the AL2023 package is older), tmux, and a data volume.
set -euxo pipefail
dnf install -y tmux unzip
cd /tmp
curl -sSfL "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" -o awscliv2.zip
unzip -q awscliv2.zip && ./aws/install --update
/usr/local/bin/aws --version
mkdir -p /data && chown ec2-user:ec2-user /data
touch /var/tmp/userdata_done
