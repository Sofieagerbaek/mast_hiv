#!/bin/bash
# Run OpenRDP's 7 methods x 3 thresholds as 21 independent single-core jobs,
# in parallel (each openrdp invocation is single-threaded; mpi4py is not
# installed so openrdp's own internal parallel path is inactive). Then merge
# each threshold's 7 per-method CSVs into one combined CSV matching the
# format the repo's run_openrdp.sh would have produced in one call.
set -u
cd "$(dirname "$0")"   # rdp_analysis/

DATA="../data"
CFG="config/hiv_rdp.ini"
METHODS=(rdp geneconv bootscan maxchi chimaera siscan threeseq)
THRESHOLDS=(0.1 0.5 0.9)
LOGDIR="results/logs"
mkdir -p "$LOGDIR"

echo "=== [$(date)] launching $(( ${#METHODS[@]} * ${#THRESHOLDS[@]} )) jobs ===" | tee -a "$LOGDIR/parallel_driver.log"

pids=()
for THR in "${THRESHOLDS[@]}"; do
  ALN="$DATA/${THR}_trimmed_alns_HIV_2021_subset.fasta"
  mkdir -p "results/${THR}/parts"
  for M in "${METHODS[@]}"; do
    OUT="results/${THR}/parts/openrdp_${THR}_${M}.csv"
    LOG="$LOGDIR/openrdp_${THR}_${M}.log"
    (
      echo "=== [$(date)] START ${THR}/${M} ==="
      openrdp "$ALN" -c "$CFG" -o "$OUT" -m "$M" -v
      echo "=== [$(date)] END ${THR}/${M} (exit $?) ==="
    ) > "$LOG" 2>&1 &
    pids+=($!)
  done
done

echo "=== [$(date)] waiting on ${#pids[@]} jobs: ${pids[*]} ===" | tee -a "$LOGDIR/parallel_driver.log"
wait
echo "=== [$(date)] all 21 jobs finished ===" | tee -a "$LOGDIR/parallel_driver.log"

# Merge per-method parts into one combined CSV per threshold.
for THR in "${THRESHOLDS[@]}"; do
  OUT="results/${THR}/openrdp_${THR}.csv"
  first=1
  : > "$OUT"
  for M in "${METHODS[@]}"; do
    PART="results/${THR}/parts/openrdp_${THR}_${M}.csv"
    if [[ ! -f "$PART" ]]; then
      echo "MISSING part: $PART" | tee -a "$LOGDIR/parallel_driver.log"
      continue
    fi
    if [[ $first -eq 1 ]]; then
      cat "$PART" >> "$OUT"
      first=0
    else
      tail -n +2 "$PART" >> "$OUT"
    fi
  done
  echo "=== [$(date)] merged -> $OUT ($(wc -l < "$OUT") lines) ===" | tee -a "$LOGDIR/parallel_driver.log"
done

echo "=== [$(date)] ALL DONE ===" | tee -a "$LOGDIR/parallel_driver.log"
