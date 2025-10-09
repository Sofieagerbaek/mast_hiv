#!/bin/bash
# Run site concordance factors (SCF) for trees with IQ-TREE
# Usage:
#   ./run_scf.sh all             # run SCF for all trees in TREES file
#   ./run_scf.sh 31 8 30         # run SCF only for trees 31, 8, and 30 

# Folders
FASTA=data/18_ML_trees/2708_trimmed_alns_HIV_2021_subset.fasta
TREES=data/35_trimmed_ML_trees.newick
SINGLETREES=data/ind_trees
OUTDIR=output/scf
SCF_REPS=100

mkdir -p "$OUTDIR"
mkdir -p "$SINGLETREES"

# To run on all 
if [ "$1" == "all" ]; then
    NUM_TREES=$(grep -c ";" "$TREES")  # count trees in file
    TREELIST=$(seq 1 $NUM_TREES) #call for each of the trees 
else
    TREELIST=$@ # otherwise tree-list is the numbers given in the commandline
fi

# Loop over selected trees
for i in $TREELIST; do
    echo "Extracting tree $i..."
    TREEFILE="$SINGLETREES/tree_${i}.newick" #putting each single tree file in the folder in data/
    sed -n "${i}p" "$TREES" > "$TREEFILE" 

    echo "Running IQ-TREE SCF for tree $i..."
    iqtree2 -te "$TREEFILE" \
            -s "$FASTA" \
            --scfl $SCF_REPS \
            --prefix "$OUTDIR/scf_tree_${i}" \
            -nt AUTO
done

echo "scf done"