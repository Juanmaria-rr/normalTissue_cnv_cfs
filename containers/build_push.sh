#!/usr/bin/env bash
# Build the CNV_cfs images for linux/amd64 and push them to ECR (eu-west-1).
# Creates each repository on first use (immutable tags, AES256, HealthOmics pull policy,
# the same settings as the account's other HealthOmics repositories).
#
# Tag = <version>-<git short sha>, so an image always points to the commit that built it.
#
# Usage: bash containers/build_push.sh [align|copynumber|all] [--no-push]
set -euo pipefail

WHICH="${1:-all}"
PUSH=1; [[ "${2:-}" == "--no-push" ]] && PUSH=0
PROFILE=<AWS_PROFILE>
REGION=eu-west-1
ACCOUNT=<AWS_ACCOUNT_ID>
VERSION=0.1.0
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SHA=$(git -C "$ROOT" rev-parse --short HEAD)
REGISTRY="$ACCOUNT.dkr.ecr.$REGION.amazonaws.com"
AWS=(aws --profile "$PROFILE" --region "$REGION")

POLICY='{"Version":"2012-10-17","Statement":[{"Sid":"AllowHealthOmicsPull","Effect":"Allow","Principal":{"Service":"omics.amazonaws.com"},"Action":["ecr:BatchGetImage","ecr:GetDownloadUrlForLayer","ecr:BatchCheckLayerAvailability"]}]}'

ensure_repo() {
  local repo="$1"
  if ! "${AWS[@]}" ecr describe-repositories --repository-names "$repo" >/dev/null 2>&1; then
    "${AWS[@]}" ecr create-repository --repository-name "$repo" \
      --image-tag-mutability IMMUTABLE --encryption-configuration encryptionType=AES256 \
      --tags Key=Project,Value=cnv_cfs >/dev/null
    "${AWS[@]}" ecr set-repository-policy --repository-name "$repo" --policy-text "$POLICY" >/dev/null
    echo "created ECR repository $repo"
  fi
}

build() {
  local name="$1"
  local repo="cnv-cfs/$name" tag="$VERSION-$SHA"
  if [[ -n "$(git -C "$ROOT" status --porcelain -- containers external)" ]]; then
    echo "uncommitted changes under containers/ or external/: commit first so the tag is truthful" >&2
    exit 1
  fi
  # No provenance/SBOM attestations: with the containerd image store they turn the push into
  # an OCI index whose attestation manifest can end up as the (immutable) tag in ECR.
  docker buildx build --platform linux/amd64 --load --provenance=false --sbom=false \
    -f "$ROOT/containers/$name/Dockerfile" -t "$REGISTRY/$repo:$tag" "$ROOT"
  echo "built $REGISTRY/$repo:$tag"
  if [[ "$PUSH" == 1 ]]; then
    ensure_repo "$repo"
    "${AWS[@]}" ecr get-login-password | docker login --username AWS --password-stdin "$REGISTRY"
    docker push "$REGISTRY/$repo:$tag"
    # The tag must resolve to the image itself, not to an attestation manifest.
    if "${AWS[@]}" ecr batch-get-image --repository-name "$repo" --image-ids imageTag="$tag" \
         --query 'images[0].imageManifest' --output text | grep -q 'in-toto'; then
      echo "ERROR: $repo:$tag points to an attestation manifest" >&2; exit 1
    fi
    echo "pushed $REGISTRY/$repo:$tag"
  fi
}

case "$WHICH" in
  align|copynumber) build "$WHICH" ;;
  all) build align; build copynumber ;;
  *) echo "usage: $0 [align|copynumber|all] [--no-push]" >&2; exit 2 ;;
esac
