#!/usr/bin/env bash
cd "$(dirname "$0")/.." || exit 1
pwd
ls
#PROJ=/Users/gzd130/Library/CloudStorage/OneDrive-RegionSjælland/Mol_phylo/revision/mast_hiv
#export PATH=$PATH:/projects/racimolab/people/fnl270/GitHub/mast_hiv/bin
FASTA=data/0.1_trimmed_alns_HIV_2021_subset.fasta
TREES=generated/0.1/35_trimmed_ML_trees.newick
#OUT1=$PROJ/output/mast_35/mast_35_fast # changed folder and file names for the plotting scripts: mast_model_35/mast_model_35 (stupid naming inconsistency)
OUT2=output/0.1/mast_model_35/mast_model_final_35 # changed folder and file names for the plotting scripts: - mast_model_35/mast_model_final_35 (stupid naming inconsistency)
mkdir -p output/0.1/mast_model_35

#iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wspmr -wslmr -redo
iqtree3 -s $FASTA -m "GTR+FO+R4+T" -te "$TREES" --prefix "$OUT2" -nt 8 -wspmr -wslmr -redo