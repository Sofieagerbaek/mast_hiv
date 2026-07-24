here::i_am("scripts/iqtree_modelselect.R")
library(here)
library(future)
library(future.apply)
plan(multisession)
setwd("/Users/gzd130/Library/CloudStorage/OneDrive-RegionSjælland/Mol_phylo/revision/mast_hiv")

# Load weights from .iqtree
# weights <- system("grep '^Tree weights' /maps/projects/racimolab/people/fnl270/GitHub/mast_hiv/output/mast_35/mast_35_fast.iqtree | sed 's/Tree weights: //; s/,//g'", intern = TRUE)
weights <- system("grep '^Tree weights' /Users/gzd130/Library/CloudStorage/OneDrive-RegionSjælland/Mol_phylo/revision/mast_hiv/output/0.1/mast_model_35/mast_model_final_35.iqtree | sed 's/Tree weights: //; s/,//g'", intern = TRUE)
weights <- as.numeric(strsplit(weights, " +")[[1]])
tree_order <- order(weights, decreasing = T)
nwk <- as.matrix(read.delim(here("generated/0.1/35_trimmed_ML_trees.newick"), header = F))

# Create subset newick files based on order
# from -fast
# for(i in 2:34) {
#   write(nwk[tree_order[1:i]], here("data/tree_subsets/",paste0("modeltrees_",i,".newick")))
# }

dir.create(
  here("generated/0.1/tree_subsets"),
  recursive = TRUE
)

for(i in 2:34) {
  write(nwk[tree_order[1:i]], here("generated/0.1/tree_subsets/",paste0("modeltrees_final_",i,".newick")))
}

FASTA <- here("data/0.1_trimmed_alns_HIV_2021_subset.fasta")
TREES_SUBSETS <- here("generated/0.1/tree_subsets")
OUT <- here("output/0.1")
commands <- character(33L)
# with -fast
# for(i in 2:34) {
#     commands[i] <- paste0("mkdir ", OUT, "/mast_model_",i, "; ", # create output subfolders
#     "export PATH=$PATH:/projects/racimolab/people/fnl270/GitHub/mast_hiv/bin; ", # export PATH
#     "iqtree2 -s ",FASTA," -m 'GTR+FO+R4+T' -te ",TREES_SUBSETS,"/modeltrees_",i,".newick --prefix ",OUT,"/mast_model_",i,"/mast_model_",i," -fast -nt 2 -wspmr -wslmr -redo") # iqtree command
# }

for(i in 2:34) {
    commands[i] <- paste0("mkdir ", OUT, "/mast_model_",i, "; ", # create output subfolders
    #"export PATH=$PATH:/projects/racimolab/people/fnl270/GitHub/mast_hiv/bin; ", # export PATH
    "iqtree3 -s ",FASTA," -m 'GTR+FO+R4+T' -te ",TREES_SUBSETS,"/modeltrees_final_",i,".newick --prefix ",OUT,"/mast_model_",i,"/mast_model_final_",i," -nt 8 -wspmr -wslmr -redo") # iqtree command
}

commands <- commands[-1]

future_lapply(commands, function(cmd) {
    system(cmd)  # Output will be printed directly to the console
}, future.seed = T)

