#!/bin/bash
# Run site concordance factors (SCF) for trees with IQ-TREE
#
# Usage (run from project root):
#   bash scripts/iqtree_scf.sh 0.5 all           # all 35 trees, 0.5 threshold
#   bash scripts/iqtree_scf.sh 0.5 30 10 21      # specific trees only
#   bash scripts/iqtree_scf.sh 0.1 28 22 27 13   # 0.1 best-model trees
#   bash scripts/iqtree_scf.sh 0.9 5 33 25 17 24 # 0.9 best-model trees
#
# Arguments:
#   $1  gap threshold (0.1 | 0.5 | 0.9)
#   $2+ tree numbers to process, or "all"
#
# Output goes to output/{output_dir}/scf/
# Output dir mapping: 0.5 -> 0.5_output, others -> same as threshold
#
# IQ-TREE binary: defaults to the checked-in ./bin/iqtree2 (Linux build, for the
# cluster). Override with IQTREE_BIN, e.g. IQTREE_BIN=iqtree3 on macOS/arm64, where
# the repo binary won't run. Flags differ between the v2 and v3 CLIs (fixed-tree flag,
# thread flag), so the binary's reported major version picks the right ones.

THRESHOLD="${1:?Usage: iqtree_scf.sh <threshold> <tree_nrs|all>}"
shift   # remaining args are tree numbers (or "all")

IQTREE_BIN="${IQTREE_BIN:-./bin/iqtree2}"
IQTREE_MAJOR=$("$IQTREE_BIN" --version 2>&1 | grep -oE "version [0-9]+" | head -1 | grep -oE "[0-9]+")
if [ "$IQTREE_MAJOR" -ge 3 ] 2>/dev/null; then
    FIXTREE_FLAGS=(--tree-fix -t)
    THREAD_FLAG=(-T AUTO)
else
    FIXTREE_FLAGS=(-te)
    THREAD_FLAG=(-nt AUTO)
fi

# Resolve output directory name (0.5 uses legacy "0.5_output" folder)
if [ "$THRESHOLD" == "0.5" ]; then
    OUTBASE="output/0.5_output"
else
    OUTBASE="output/${THRESHOLD}"
fi

FASTA="data/${THRESHOLD}_trimmed_alns_HIV_2021_subset.fasta"
TREES="generated/${THRESHOLD}/35_trimmed_ML_trees.newick"
SINGLETREES="generated/${THRESHOLD}/ind_trees"
OUTDIR="${OUTBASE}/scf"
SCF_REPS=100

# Validate inputs
if [ ! -f "$FASTA" ]; then
    echo "ERROR: alignment not found: $FASTA" >&2; exit 1
fi
if [ ! -f "$TREES" ]; then
    echo "ERROR: tree file not found: $TREES" >&2; exit 1
fi

mkdir -p "$OUTDIR"
mkdir -p "$SINGLETREES"

# Build list of tree indices to process
if [ "$1" == "all" ]; then
    NUM_TREES=$(grep -c ";" "$TREES")
    TREELIST=$(seq 1 "$NUM_TREES")
else
    TREELIST="$@"
fi

# Loop over selected trees
for i in $TREELIST; do
    echo "--- Tree ${i} (threshold ${THRESHOLD}) ---"
    TREEFILE="${SINGLETREES}/tree_${i}.newick"

    echo "  Extracting tree ${i}..."
    sed -n "${i}p" "$TREES" > "$TREEFILE"

    echo "  Running IQ-TREE SCF..."
    "$IQTREE_BIN" \
        "${FIXTREE_FLAGS[@]}" "$TREEFILE" \
        -s   "$FASTA" \
        --scfl "$SCF_REPS" \
        --prefix "${OUTDIR}/scf_tree_${i}" \
        "${THREAD_FLAG[@]}"
done

echo "SCF done for threshold ${THRESHOLD}."
