#!/bin/bash
PROJ=/projects/racimolab/people/fnl270/GitHub/mast_hiv
DATA=$PROJ/data
export PATH=$PATH:/projects/racimolab/people/fnl270/GitHub/mast_hiv/bin
FASTA=$DATA/18_ML_trees/2708_trimmed_alns_HIV_2021_subset.fasta
TREES=$DATA/35_trimmed_ML_trees.newick
OUT1=$PROJ/output/mast_35/mast_35_fast # changed folder and file names for the plotting scripts: mast_model_35/mast_model_35 (stupid naming inconsistency)
OUT2=$PROJ/output/mast_35/mast_35 # changed folder and file names for the plotting scripts: - mast_model_35/mast_model_final_35 (stupid naming inconsistency)

iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wspmr -wslmr -redo
iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT2 -nt AUTO -wspmr -wslmr -redo