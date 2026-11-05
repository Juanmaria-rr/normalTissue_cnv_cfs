#!/usr/bin/env bash
# Copy SRA run files from the public SRA open-data bucket (us-east-1, read unsigned)
# to our eu-west-1 bucket, where HealthOmics can read them. Runs on the staging EC2.
#
# Each run is downloaded to --workdir, checked against the source size, then the
# whole directory is synced to --dest with the signed profile. Re-running skips runs
# already downloaded with the right size. Writes a manifest (run, bytes, sha256) that is
# uploaded next to the data.
#
# Usage:
#   bash stage_sra_to_s3.sh --runs runs_pilot.txt \
#        --dest s3://<PROJECT_BUCKET>/cnv_cfs/sra \
#        --profile <AWS_PROFILE> --workdir /data/sra [--jobs 32]
set -euo pipefail

SRC_BUCKET="s3://sra-pub-run-odp/sra"
JOBS=32
while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs)    RUNS="$2"; shift 2 ;;
    --dest)    DEST="${2%/}"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --workdir) WORKDIR="${2%/}"; shift 2 ;;
    --jobs)    JOBS="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
: "${RUNS:?--runs required}" "${DEST:?--dest required}" "${PROFILE:?--profile required}" "${WORKDIR:?--workdir required}"

mkdir -p "$WORKDIR" "$WORKDIR/../logs"
LOG="$WORKDIR/../logs/stage_$(date -u +%Y%m%dT%H%M%SZ)"
aws sts get-caller-identity --profile "$PROFILE" --query Arn --output text   # fail fast

# --- 1. download (unsigned), one run per job --------------------------------------
fetch_one() {
  local run="$1" workdir="$2"
  local out="$workdir/$run/$run"
  local expected
  expected=$(aws s3api head-object --no-sign-request --region us-east-1 \
             --bucket sra-pub-run-odp --key "sra/$run/$run" --query ContentLength --output text \
             2>/dev/null) || true
  if [[ -z "$expected" || "$expected" == "None" ]]; then
    echo -e "$run\tMISSING_AT_SOURCE"; return 0
  fi
  if [[ -f "$out" && $(stat -c %s "$out") == "$expected" ]]; then
    echo -e "$run\tSKIP\t$expected"; return 0
  fi
  mkdir -p "$workdir/$run"
  if aws s3 cp --no-sign-request --region us-east-1 --only-show-errors \
       "$SRC_BUCKET/$run/$run" "$out"; then
    local got; got=$(stat -c %s "$out")
    if [[ "$got" == "$expected" ]]; then echo -e "$run\tOK\t$got"
    else echo -e "$run\tSIZE_MISMATCH\t$got\t$expected"; fi
  else
    echo -e "$run\tDOWNLOAD_FAILED"
  fi
}
export -f fetch_one
export SRC_BUCKET

grep -v '^\s*$' "$RUNS" | sort -u > "$LOG.runs"
N=$(wc -l < "$LOG.runs")
echo "[$(date -u +%T)] downloading $N runs with $JOBS jobs"
xargs -P "$JOBS" -I{} bash -c 'fetch_one "$@"' _ {} "$WORKDIR" \
  < "$LOG.runs" > "$LOG.download.tsv"
awk -F'\t' '{print $2}' "$LOG.download.tsv" | sort | uniq -c
BAD=$(awk -F'\t' '$2!="OK" && $2!="SKIP"' "$LOG.download.tsv" | wc -l)
if [[ "$BAD" -gt 0 ]]; then
  echo "WARNING: $BAD runs not downloaded correctly (see $LOG.download.tsv); they will not be uploaded"
fi

# --- 2. manifest ----------------------------------------------------------------------
echo "[$(date -u +%T)] computing sha256"
awk -F'\t' '$2=="OK"||$2=="SKIP"{print $1}' "$LOG.download.tsv" |
  xargs -P "$JOBS" -I{} sh -c 'f="$1/$2/$2"; printf "%s\t%s\t%s\n" "$2" "$(stat -c %s "$f")" "$(sha256sum "$f" | cut -d" " -f1)"' _ "$WORKDIR" {} |
  sort > "$LOG.manifest.tsv"
sed -i '1i run\tbytes\tsha256' "$LOG.manifest.tsv"

# --- 3. upload (signed) ---------------------------------------------------------------
echo "[$(date -u +%T)] uploading to $DEST"
aws s3 sync "$WORKDIR" "$DEST" --profile "$PROFILE" --region eu-west-1 --only-show-errors
aws s3 cp "$LOG.manifest.tsv" "$DEST/_manifest/$(basename "$LOG").manifest.tsv" \
  --profile "$PROFILE" --region eu-west-1 --only-show-errors

# --- 4. verify: every manifest run present at destination with the same size ---------
echo "[$(date -u +%T)] verifying destination"
aws s3 ls "$DEST/" --recursive --profile "$PROFILE" --region eu-west-1 |
  awk '{n=split($4,p,"/"); if (p[n]==p[n-1]) print p[n]"\t"$3}' | sort > "$LOG.dest.tsv"
# run<TAB>bytes pairs in the manifest but not at the destination = missing or wrong size
MISSING=$(comm -23 <(tail -n +2 "$LOG.manifest.tsv" | cut -f1,2 | sort) "$LOG.dest.tsv" | wc -l)
echo "runs in manifest: $(($(wc -l < "$LOG.manifest.tsv") - 1))  objects at dest: $(wc -l < "$LOG.dest.tsv")  missing/size-mismatch: $MISSING"
[[ "$MISSING" -eq 0 ]] && echo "STAGING VERIFIED" || { echo "STAGING INCOMPLETE"; exit 1; }
