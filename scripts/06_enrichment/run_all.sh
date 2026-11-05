#!/usr/bin/env bash
# Rebuild the pilot enrichment analysis end to end (depth recalibration, event tables, tests A-C,
# sensitivity and robustness runs, power projection, master tables, funnel, figures, results
# document). Run from the repository root.
set -euo pipefail
R=results/enrichment/pilot
CFS=data/features/cfs/cfs_catalogs_hg19.tsv
GENES=data/features/genes
RT=external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv
TA="python3 scripts/06_enrichment/test_a_permutation.py --cfs $CFS --genes $GENES"
python3 scripts/03_cn_qc/recalibrate_depth.py --qc results/pilot_cn/qc
python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc --profiles results/pilot_cn/review --rt $RT --out $R
# robustness event sets: linkage distance and collapsing on E2-E4 breakends only
python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc --profiles results/pilot_cn/review --rt $RT --out $R --sets main --link-bp 1000000 --tag _link1mb
python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc --profiles results/pilot_cn/review --rt $RT --out $R --sets main --link-bp 5000000 --tag _link5mb
python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc --profiles results/pilot_cn/review --rt $RT --out $R --sets main --collapse-on main --tag _collapse_main
$TA --events $R/main --out $R/main
for s in nocuts samecells cells500 main_link1mb main_link5mb main_collapse_main; do
  $TA --events $R/$s --out $R/$s --perm 0
done
$TA --events $R/refined --out $R/refined --bins 200 --perm 0
for d in P07 P08 P18 P20 P36; do
  $TA --events $R/main --out $R/main --perm 0 --exclude-donor $d --tag _without_$d
done
python3 scripts/06_enrichment/test_b_regression.py --events $R/main --rt $RT --cfs $CFS --genes $GENES --out $R/main
python3 scripts/06_enrichment/test_c_per_donor.py --events $R/main --cfs $CFS --genes $GENES --out $R/main
python3 scripts/06_enrichment/robustness.py --res $R
python3 scripts/06_enrichment/power_scaleup.py --res $R --out $R/power_scaleup.tsv
python3 scripts/06_enrichment/build_master_tables.py --res $R --qc results/pilot_cn/qc --align data/pilot/align_qc_pilot.tsv --cells data/metadata/cells.tsv --cfs $CFS --genes $GENES
python3 scripts/06_enrichment/funnel.py --res $R --cells data/metadata/cells.tsv --align data/pilot/align_qc_pilot.tsv --fig-dir docs/enrichment
python3 scripts/docs/enrichment_figures.py --res $R --out docs/enrichment
python3 scripts/docs/build_enrichment_html.py --res $R --out docs/enrichment
python3 scripts/docs/build_tests_explainer.py --res $R --out docs/enrichment
python3 scripts/docs/build_enrichment_brief.py --res $R --out docs/enrichment --standalone
echo "all done"
