#!/bin/zsh
# Run the MAST tree-number power bootstrap for all three gap thresholds.
# Run from the repo root (here::i_am anchors on mast_hiv.Rproj):
#     zsh power_bootstrap/run_power_bootstrap.sh
#
# Override scale for a quick check, e.g.:
#     PB_B=2 PB_KMAX=4 zsh power_bootstrap/run_power_bootstrap.sh
set -u
LOGDIR="power_bootstrap/results/logs"
mkdir -p "$LOGDIR"

for THR in 0.1 0.5 0.9; do
  echo "=== [$(date)] START threshold $THR ==="
  Rscript power_bootstrap/power_bootstrap.R "$THR" >"$LOGDIR/power_${THR}.log" 2>&1
  echo "=== [$(date)] END threshold $THR (exit $?) ==="
done
echo "=== [$(date)] ALL DONE ==="
