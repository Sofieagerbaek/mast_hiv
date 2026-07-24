#!/bin/zsh
# Run OpenRDP on the three gap-trimmed HIV subtype consensus alignments.
#
# OpenRDP is NOT installed in this environment yet. Install it first, for
# example with:
#     python3 -m pip install openrdp
# or from a clone of https://github.com/PoonLab/OpenRDP:
#     python3 -m pip install --user .
#
# Then run this script from the rdp_analysis folder:
#     cd mast_hiv/rdp_analysis && zsh scripts/run_openrdp.sh
#
# Each alignment is scanned with all seven methods and its events are written
# to results/<threshold>/openrdp_<threshold>.csv.

set -u

HERE="${0:A:h}/.."            # rdp_analysis/
DATA="$HERE/../data"          # mast_hiv/data
CFG="$HERE/config/hiv_rdp.ini"
METHODS=(rdp geneconv bootscan maxchi chimaera siscan threeseq)

if ! command -v openrdp >/dev/null 2>&1; then
  echo "ERROR: 'openrdp' not found on PATH. Install it first (see header)." >&2
  exit 1
fi

for THR in 0.1 0.5 0.9; do
  ALN="$DATA/${THR}_trimmed_alns_HIV_2021_subset.fasta"
  OUTDIR="$HERE/results/${THR}"
  OUT="$OUTDIR/openrdp_${THR}.csv"
  mkdir -p "$OUTDIR"

  if [[ ! -s "$ALN" ]]; then
    echo "SKIP ${THR}: alignment not found at $ALN" >&2
    continue
  fi

  echo "=== [$(date)] OpenRDP on ${THR} -> ${OUT} ==="
  openrdp "$ALN" -c "$CFG" -o "$OUT" -m ${METHODS[@]}
  echo "=== [$(date)] done ${THR} (exit $?) ==="
done
