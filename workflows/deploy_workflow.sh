#!/usr/bin/env bash
# Package a workflow directory under workflows/ and register it as a private HealthOmics
# workflow (eu-west-1). The workflow name carries the git short sha of the commit it was
# built from, so a run always points to traceable code.
#
# Usage: bash workflows/deploy_workflow.sh <workflow_dir_name>
#   e.g. bash workflows/deploy_workflow.sh build_reference
# Prints the new workflow id.
set -euo pipefail

NAME="$1"
PROFILE=<AWS_PROFILE>
REGION=eu-west-1
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/workflows/$NAME"
[[ -f "$DIR/main.nf" ]] || { echo "no $DIR/main.nf" >&2; exit 2; }
if [[ -n "$(git -C "$ROOT" status --porcelain -- "workflows/$NAME")" ]]; then
  echo "uncommitted changes in workflows/$NAME: commit first" >&2; exit 1
fi
SHA=$(git -C "$ROOT" rev-parse --short HEAD)
ZIP=$(mktemp -d)/"$NAME.zip"
(cd "$DIR" && zip -qr "$ZIP" . -x 'parameter-template.json')

aws omics create-workflow --profile "$PROFILE" --region "$REGION" \
  --name "cnv-cfs-$NAME-$SHA" --engine NEXTFLOW --main main.nf \
  --definition-zip "fileb://$ZIP" \
  --parameter-template "file://$DIR/parameter-template.json" \
  --storage-type DYNAMIC \
  --tags Project=cnv_cfs,GitSha="$SHA" \
  --query id --output text
