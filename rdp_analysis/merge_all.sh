#!/bin/bash
# Rebuild results/<threshold>/openrdp_<threshold>.csv from whatever per-method
# files are present in results/<threshold>/parts/.
#
# Run this after rerun_siscan.sh finishes, to fold SiScan into the combined CSVs.
# Safe to run at any time: it reports which methods contributed and which are
# missing, so a partial analysis is visible rather than looking complete.
set -u
cd "$(dirname "$0")"   # rdp_analysis/

METHODS=(rdp geneconv bootscan maxchi chimaera siscan threeseq)
HEADER='Method,Start,End,Recombinant,Parent1,Parent2,Pvalue'

for THR in 0.1 0.5 0.9; do
  OUT="results/${THR}/openrdp_${THR}.csv"
  TMP="${OUT}.tmp"
  echo "$HEADER" > "$TMP"
  present=() missing=()
  for M in "${METHODS[@]}"; do
    PART="results/${THR}/parts/openrdp_${THR}_${M}.csv"
    if [[ ! -s "$PART" ]]; then
      missing+=("$M")
      continue
    fi
    n=$(( $(wc -l < "$PART") - 1 ))
    tail -n +2 "$PART" >> "$TMP"
    present+=("${M}:${n}")
  done
  mv "$TMP" "$OUT"
  echo "${THR}: $(( $(wc -l < "$OUT") - 1 )) rows  [${present[*]}]"
  if ((${#missing[@]})); then
    echo "      no output from: ${missing[*]}"
  fi
done
