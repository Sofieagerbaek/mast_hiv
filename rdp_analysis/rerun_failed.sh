#!/bin/bash
# Re-run the OpenRDP methods that produced no usable output in the first pass
# (see results/logs/parallel_driver.log, 2026-07-24 -> 2026-07-27):
#
#   maxchi    crashed on all three thresholds inside openrdp itself
#   chimaera  ran clean but reported zero events on all three thresholds
#   siscan    killed by the kernel partway through 0.1 and 0.5
#
# maxchi and chimaera are re-run because of the local openrdp patches in
# patches/ (a wrong validity guard in common.calculate_chi2 that both crashed
# maxchi and made both methods silently discard perfectly-associated windows,
# plus an off-by-one when reading the chi2 array). siscan is unaffected by those
# patches -- it never calls calculate_chi2 -- so its 0.9 result (zero events)
# stands and only 0.1 and 0.5 are re-run, this time split across processes by
# siscan_chunk.py so each one stays small and finishes in hours not days.
#
# Unlike the first-pass script this one captures exit statuses correctly. The
# original wrote `echo "... ($(date)) ... (exit $?)"`, and the $(date) command
# substitution resets $? before it is expanded, so every job was logged as
# "exit 0" including the ones that crashed or were killed.
set -u
cd "$(dirname "$0")"   # rdp_analysis/

# Refuse to run against an unpatched openrdp: unpatched, MaxChi crashes, SiScan
# loops forever leaking ~5 GB/h, and MaxChi/Chimaera silently disagree with the
# committed results. See check_openrdp_patches.py.
python3 check_openrdp_patches.py || exit 1


DATA="../data"
CFG="config/hiv_rdp.ini"
LOGDIR="results/logs/rerun"
DRIVER="$LOGDIR/rerun_driver.log"
NCHUNKS="${SISCAN_NCHUNKS:-8}"
mkdir -p "$LOGDIR"

log() { echo "=== [$(date)] $* ===" | tee -a "$DRIVER"; }

# Run a command, log its real exit status, keep the log path stable.
run_job() {
  local name="$1"; shift
  local logf="$LOGDIR/${name}.log"
  {
    echo "=== START ${name} $(date) ==="
    "$@"
    rc=$?          # captured before any command substitution can clobber it
    echo "=== END ${name} $(date) (exit ${rc}) ==="
    exit $rc
  } > "$logf" 2>&1
}

log "rerun starting (siscan split into ${NCHUNKS} chunks per threshold)"

declare -A JOBNAME
pids=()

# ---- maxchi and chimaera: whole-alignment, one core each -------------------
for THR in 0.1 0.5 0.9; do
  ALN="$DATA/${THR}_trimmed_alns_HIV_2021_subset.fasta"
  for M in maxchi chimaera; do
    mkdir -p "results/${THR}/parts"
    OUT="results/${THR}/parts/openrdp_${THR}_${M}.csv"
    run_job "openrdp_${THR}_${M}" \
      openrdp "$ALN" -c "$CFG" -o "$OUT" -m "$M" -v &
    JOBNAME[$!]="${THR}/${M}"
    pids+=($!)
  done
done

# ---- siscan: 0.1 and 0.5 only, chunked ------------------------------------
for THR in 0.1 0.5; do
  ALN="$DATA/${THR}_trimmed_alns_HIV_2021_subset.fasta"
  mkdir -p "results/${THR}/parts/siscan_chunks"
  for ((c = 0; c < NCHUNKS; c++)); do
    OUT="results/${THR}/parts/siscan_chunks/siscan_${THR}_chunk${c}.pkl"
    run_job "openrdp_${THR}_siscan_chunk${c}" \
      python3 siscan_chunk.py "$ALN" -c "$CFG" -o "$OUT" \
        --chunk "$c" --nchunks "$NCHUNKS" &
    JOBNAME[$!]="${THR}/siscan[${c}]"
    pids+=($!)
  done
done

log "launched ${#pids[@]} jobs: ${pids[*]}"

# ---- memory sampler ------------------------------------------------------
# The first-pass siscan jobs were killed with no record of why (dmesg needs root
# here). Sample RSS every 5 min so a second kill leaves a usable trace.
(
  while :; do
    ts=$(date '+%Y-%m-%dT%H:%M:%S')
    tot=$(ps -o rss= -C python3 -C openrdp 2>/dev/null | awk '{s+=$1} END {printf "%.1f", s/1048576}')
    avail=$(awk '/MemAvailable/ {printf "%.1f", $2/1048576}' /proc/meminfo)
    top=$(ps -o rss=,cmd= -C python3 2>/dev/null | sort -rn | head -1 | awk '{printf "%.2fGB", $1/1048576}')
    echo "$ts rss_total=${tot}GB mem_available=${avail}GB largest_python=${top}"
    sleep 300
  done
) > "$LOGDIR/memory_sampler.log" 2>&1 &
SAMPLER=$!

# ---- wait, reporting each job's real status ------------------------------
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

log "all jobs finished (${failed} failed)"

# ---- merge siscan chunks -------------------------------------------------
for THR in 0.1 0.5; do
  CHUNKDIR="results/${THR}/parts/siscan_chunks"
  if python3 merge_siscan_chunks.py "$CHUNKDIR"/siscan_"${THR}"_chunk*.pkl \
       --nchunks "$NCHUNKS" \
       -o "results/${THR}/parts/openrdp_${THR}_siscan.csv" >> "$DRIVER" 2>&1; then
    log "merged siscan chunks -> results/${THR}/parts/openrdp_${THR}_siscan.csv"
  else
    log "FAIL merging siscan chunks for ${THR} (see $DRIVER)"
    failed=$((failed + 1))
  fi
done

# ---- re-merge the combined per-threshold CSVs ----------------------------
METHODS=(rdp geneconv bootscan maxchi chimaera siscan threeseq)
for THR in 0.1 0.5 0.9; do
  OUT="results/${THR}/openrdp_${THR}.csv"
  echo 'Method,Start,End,Recombinant,Parent1,Parent2,Pvalue' > "$OUT"
  for M in "${METHODS[@]}"; do
    PART="results/${THR}/parts/openrdp_${THR}_${M}.csv"
    if [[ ! -s "$PART" ]]; then
      log "NOTE ${THR}/${M}: no part file, contributing 0 rows"
      continue
    fi
    tail -n +2 "$PART" >> "$OUT"
  done
  log "merged -> $OUT ($(( $(wc -l < "$OUT") - 1 )) data rows)"
done

log "RERUN DONE (${failed} failures)"
