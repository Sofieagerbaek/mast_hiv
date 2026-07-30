#!/bin/bash
# Run the MAST tree-number power bootstrap for all three gap thresholds,
# sequentially. Run from the repo root (here::i_am anchors on mast_hiv.Rproj):
#
#     bash power_bootstrap/run_power_bootstrap.sh
#
# The thresholds run one after another rather than all at once, so the core cap
# below applies to one threshold at a time and the machine stays usable. Each
# threshold's own k = 2..10 refits are what run in parallel, PB_NCORES at a time.
#
# Point it at your IQ-TREE 3 and cap the cores on a shared machine:
#
#     PB_IQTREE=/path/to/iqtree3 PB_NCORES=20 bash power_bootstrap/run_power_bootstrap.sh
#
# Override scale for a quick end-to-end check:
#
#     PB_B=2 PB_KMAX=4 PB_NCORES=4 bash power_bootstrap/run_power_bootstrap.sh
#
# The committed results in results/<t>/ were produced with:
#
#     PB_IQTREE=~/iqtree3 PB_NCORES=20 bash power_bootstrap/run_power_bootstrap.sh
#
# at PB_B=100, PB_KMIN=2, PB_KMAX=10, PB_SEED=1 (all defaults), on a 64-core
# Linux server, taking ~5.2 days of wall clock in total.
set -u

LOGDIR="power_bootstrap/results/logs"
DRIVER_LOG="$LOGDIR/sequential_driver.log"
mkdir -p "$LOGDIR"

for THR in 0.1 0.5 0.9; do
  echo "=== [$(date)] START threshold $THR ===" | tee -a "$DRIVER_LOG"
  Rscript power_bootstrap/power_bootstrap.R "$THR" >"$LOGDIR/power_${THR}.log" 2>&1
  # Capture the status before anything else can touch $?. An earlier version of
  # this script wrote `echo "... $(date) ... (exit $?)"`, where the command
  # substitution resets $? before it is expanded, so every threshold logged as
  # "exit 0" whether it succeeded or not.
  STATUS=$?
  echo "=== [$(date)] END threshold $THR (exit $STATUS) ===" | tee -a "$DRIVER_LOG"
done
echo "=== [$(date)] ALL DONE ===" | tee -a "$DRIVER_LOG"
