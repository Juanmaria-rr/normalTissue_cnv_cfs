#!/usr/bin/env bash
# Launch the curation app against one batch.
#
#   ./run_curation.sh              interactive menu
#   ./run_curation.sh batch_01     straight to a batch
#
# Curate ONE batch per run. The app scans for samples at startup, so stop it
# with Ctrl-C and relaunch to switch batches.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$HERE/.venv"
PORT="${CURATION_PORT:-8501}"

if [ ! -x "$VENV/bin/shiny" ]; then
    echo "ERROR: virtualenv not found. Run ./setup_mac.sh first." >&2
    exit 1
fi

# Batch list comes from the manifest, so the launcher never hardcodes labels.
LABELS=()
FLAGS=()
while IFS=$'\t' read -r lab flag || [ -n "$lab" ]; do
    [ -z "$lab" ] && continue
    LABELS+=("$lab")
    FLAGS+=("$flag")
done < "$HERE/data/batches.tsv"

count_samples() {
    local f="$HERE/data/$1/samples.txt"
    if [ -f "$f" ]; then grep -c . "$f"; else echo "?"; fi
}

choose() {
    echo "Curation - pick a batch"
    echo "======================="
    echo
    local i=1
    for idx in "${!LABELS[@]}"; do
        local mark=""
        [ -n "${FLAGS[$idx]}" ] && mark="  <- ${FLAGS[$idx]}"
        printf "  %d) %-10s %3s samples%s\n" "$i" "${LABELS[$idx]}" \
               "$(count_samples "${LABELS[$idx]}")" "$mark"
        i=$((i + 1))
    done
    echo
    read -r -p "Choice [1-${#LABELS[@]}]: " n
    if ! [[ "$n" =~ ^[0-9]+$ ]] || [ "$n" -lt 1 ] || [ "$n" -gt "${#LABELS[@]}" ]; then
        echo "Not a valid choice." >&2
        exit 1
    fi
    BATCH="${LABELS[$((n - 1))]}"
}

if [ $# -ge 1 ]; then BATCH="$1"; else choose; fi

DATA="$HERE/data/$BATCH"
if [ ! -d "$DATA" ]; then
    echo "ERROR: no such batch '$BATCH'." >&2
    echo "Available: ${LABELS[*]}" >&2
    exit 1
fi

FITS="$DATA/curated_fits.csv"
echo
echo "Batch          : $BATCH ($(count_samples "$BATCH") samples)"
echo "Saves fits to  : $FITS"
if [ -f "$FITS" ]; then
    echo "                 (already holds $(( $(grep -c . "$FITS") - 1 )) curated samples)"
fi
echo "Open           : http://127.0.0.1:$PORT"
echo
echo "Press Ctrl-C to stop."
echo

cd "$HERE/app"
export CURATION_DATA_PATH="$DATA"
exec "$VENV/bin/shiny" run curation_app/app.py --port "$PORT"
