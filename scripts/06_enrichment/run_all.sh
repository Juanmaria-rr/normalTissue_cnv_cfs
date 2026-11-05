#!/usr/bin/env bash
# Rebuild the pilot enrichment analysis end to end (event tables, tests A-C, sensitivity runs,
# master tables, funnel, figures, results document). Run from the repository root.
set -euo pipefail
R=results/enrichment/pilot
CFS=data/features/cfs/cfs_catalogs_hg19.tsv
GENES=data/features/genes
RT=external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv
python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc --profiles results/pilot_cn/review --rt $RT --out $R
python3 scripts/06_enrichment/test_a_permutation.py --events $R/main --cfs $CFS --genes $GENES --out $R/main --shuffle all
python3 scripts/06_enrichment/test_a_permutation.py --events $R/nocuts --cfs $CFS --genes $GENES --out $R/nocuts --shuffle none
python3 scripts/06_enrichment/test_a_permutation.py --events $R/samecells --cfs $CFS --genes $GENES --out $R/samecells --shuffle none
for d in P07 P08 P18 P20 P36; do
  python3 scripts/06_enrichment/test_a_permutation.py --events $R/main --cfs $CFS --genes $GENES --out $R/main --shuffle none --exclude-donor $d --tag _without_$d
done
python3 scripts/06_enrichment/test_b_regression.py --events $R/main --rt $RT --cfs $CFS --genes $GENES --out $R/main
python3 scripts/06_enrichment/test_c_per_donor.py --events $R/main --cfs $CFS --genes $GENES --out $R/main
python3 scripts/06_enrichment/build_master_tables.py --res $R --qc results/pilot_cn/qc --align data/pilot/align_qc_pilot.tsv --cells data/metadata/cells.tsv --cfs $CFS --genes $GENES
python3 scripts/06_enrichment/funnel.py --res $R --cells data/metadata/cells.tsv --align data/pilot/align_qc_pilot.tsv --fig-dir docs/enrichment
python3 scripts/docs/enrichment_figures.py --res $R --out docs/enrichment
python3 scripts/docs/build_enrichment_html.py --res $R --out docs/enrichment
echo "all done"
