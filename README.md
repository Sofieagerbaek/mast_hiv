# Examining recombination in HIV M-group subtypes with IQ-tree's mixture across sites and trees (MAST) model

Scripts and files used for the analysis of a subset of 37 of the 2021 HIV subtype 1 consensus sequences (https://www.hiv.lanl.gov/content/index) using mixture models. 

- Scripts in the Scripts folder
- Output from the full MAST model run on all 35 initial trees and MAST output from the optimal 6-tree model in the ouput folder

#### Analysis steps: 

Trimming of the alignment and tree building was done with the "make_ML_trees.Rmd" script, the 18 initial trees and 17 shifted trees are in the data/ folder. 

IQ-tree's MAST (http://www.iqtree.org/doc/Complex-Models) model was run on the 36 newick trees in "35_trimmed_ML_trees.newick" and the trimmed alignment "2708_trimmed_alns_HIV_2021_subset.fasta" in the command line with the call: 

    iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wspmr -wslmr -redo

The MAST analysis was rerun on sets of trees from k=2 to k=34 iteratively adding trees in order of their weight. Bayesian Information Criterion (BIC) from these were used to decide the optimal number of trees (=6). The MAST output files from this model run is in this repo folder 'output/mast_model_6' with the prefix 'mast_model_final_6'.

Post processing and most visualisations are performed in the Rmd  "MAST_post_processing.Rmd". The knitted HTML is also provided. 

Some other smaller visualisation scripts are present in the scripts folder, such as IC and weight plotting - as well as the shell script used to run iqtree MAST runs (first and second) and scf calculations.
