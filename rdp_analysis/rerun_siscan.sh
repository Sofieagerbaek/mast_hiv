#!/bin/bash
# Re-run SiScan on the 0.1 and 0.5 alignments, chunked across processes.
#
# Second attempt. The first attempt (results/logs/rerun/, launched 11:52) had to
# be abandoned 2 h in: 3 of its 18 chunks stopped advancing and grew by about
# 5 GB/h each. That was the same fault that got the original single-process 0.1
# and 0.5 runs OOM-killed after 23 h and 40 h -- an infinite loop in
# Siscan.find_signal, now fixed in patches/. This run uses the fixed code, so
# every chunk should finish in about 10 h rather than one of them eating the
# machine.
#
# 0.9 is not re-run: SiScan does not go through common.calculate_chi2, so none
# of the other patches touch it, and its first-pass run completed normally with
# zero events.
set -u
cd "$(dirname "$0")"   # rdp_analysis/

DATA="../data"
CFG="config/hiv_rdp.ini"
LOGDIR="results/logs/siscan_rerun2"
DRIVER="$LOGDIR/driver.log"
NCHUNKS="${SISCAN_NCHUNKS:-9}"
mkdir -p "$LOGDIR"

log() { echo "=== [$(date)] $* ===" | tee -a "$DRIVER"; }

log "starting SiScan rerun 2: thresholds 0.1 and 0.5, ${NCHUNKS} chunks each"

declare -A JOBNAME
pids=()
for THR in 0.1 0.5; do
  ALN="$DATA/${THR}_trimmed_alns_HIV_2021_subset.fasta"
  CHUNKDIR="results/${THR}/parts/siscan_chunks"
  mkdir -p "$CHUNKDIR"
  # Drop any pickles from the abandoned first attempt so a stale file cannot be
  # mistaken for a completed chunk by the merge step.
  rm -f "$CHUNKDIR"/siscan_"${THR}"_chunk*.pkl
  for ((c = 0; c < NCHUNKS; c++)); do
    logf="$LOGDIR/siscan_${THR}_chunk${c}.log"
    (
      echo "=== START ${THR}/siscan[${c}] $(date) ==="
      python3 siscan_chunk.py "$ALN" -c "$CFG" \
        -o "$CHUNKDIR/siscan_${THR}_chunk${c}.pkl" \
        --chunk "$c" --nchunks "$NCHUNKS"
      rc=$?
      echo "=== END ${THR}/siscan[${c}] $(date) (exit ${rc}) ==="
      exit $rc
    ) > "$logf" 2>&1 &
    JOBNAME[$!]="${THR}/siscan[${c}]"
    pids+=($!)
  done
done

log "launched ${#pids[@]} chunks: ${pids[*]}"

# Watch for the leak signature returning. A healthy chunk sits at ~120 MB; the
# runaway ones passed 6 GB. Flag anything over 2 GB rather than waiting for the
# OOM killer to make the decision.
(
  while :; do
    ts=$(date '+%Y-%m-%dT%H:%M:%S')
    big=$(ps -o rss=,args= -C python3 2>/dev/null | grep siscan_chunk \
          | awk '$1 > 2097152 {printf "%.1fGB:%s ", $1/1048576, $(NF-2)}')
    tot=$(ps -o rss= -C python3 2>/dev/null | awk '{s+=$1} END {printf "%.1f", s/1048576}')
    avail=$(awk '/MemAvailable/ {printf "%.1f", $2/1048576}' /proc/meminfo)
    echo "$ts rss_total=${tot}GB mem_available=${avail}GB over_2GB=[${big:-none}]"
    sleep 300
  done
) > "$LOGDIR/memory_sampler.log" 2>&1 &
SAMPLER=$!

failed=0
for pid in "${pids[@]}"; do
  if wait "$pid"; then
    log "OK   ${JOBNAME[$pid]}"
  else
    rc=$?
    log "FAIL ${JOBNAME[$pid]} (exit ${rc})"
    failed=$((failed + 1))
  fi
done
kill "$SAMPLER" 2>/dev/null

log "all chunks finished (${failed} failed)"

for THR in 0.1 0.5; do
  CHUNKDIR="results/${THR}/parts/siscan_chunks"
  if python3 merge_siscan_chunks.py "$CHUNKDIR"/siscan_"${THR}"_chunk*.pkl \
       --nchunks "$NCHUNKS" \
       -o "results/${THR}/parts/openrdp_${THR}_siscan.csv" >> "$DRIVER" 2>&1; then
    log "merged -> results/${THR}/parts/openrdp_${THR}_siscan.csv"
  else
    log "FAIL merging ${THR} (see $DRIVER)"
    failed=$((failed + 1))
  fi
done

log "SISCAN RERUN 2 DONE (${failed} failures) -- run merge_all.sh to rebuild the combined CSVs"
