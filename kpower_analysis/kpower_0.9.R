# kpower power assessment - 0.9 gap threshold dataset
#
# Assesses whether "data/0.9_trimmed_alns_HIV_2021_subset.fasta" (columns with
# gaps in >90% of sequences removed) has enough signal to reliably recover the
# best-fitting number of MAST tree classes (K), following the same k=2..16
# range explored in scripts/make_ML_trees.Rmd / iqtree_modelselect.R (optimal
# K was 6 for the full 35-tree MAST run).
#
# kpower's "+T" pathway fits its own candidate trees from alignment windows
# internally (it does not reuse the project's precomputed data/17_ML_trees
# windows), so this is a self-contained check of the alignment's power, run
# independently of the main MAST pipeline.

here::i_am("kpower_analysis/kpower_0.9.R")
library(here)
library(kpower)

THRESHOLD <- "0.9"

FASTA  <- here("data", paste0(THRESHOLD, "_trimmed_alns_HIV_2021_subset.fasta"))
OUTDIR <- here("kpower_analysis", "output", THRESHOLD)
FIGDIR <- here("kpower_analysis", "figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

check_iqtree()

# B = 1000 (the package default) is the publication-quality setting, but at
# K_max = 16 that means up to 16 * 1000 MAST refits on top of the 16 empirical
# fits, which can run for a very long time. Start small to sanity-check the
# pipeline, then raise B (and n_cores) for the real run.
result <- kpower(
  alignment  = FASTA,
  K_max      = 16,
  K_min      = 1,
  base_model = "GTR",
  mix_type   = "+T",     # MAST tree-mixture pathway
  ic         = "BIC",
  # fixed_tree = "NJ" (the kpower default) drives IQ-TREE's -t BIONJ
  # --tree-fix for the K=1 single-tree fit, which crashes
  # (computeDistanceMatrix / SIGABRT) on this machine's IQ-TREE3 3.0.1
  # (Homebrew, macOS ARM64) build regardless of thread count. Using NULL
  # falls back to a --fast heuristic search instead, which works fine.
  fixed_tree = NULL,
  B          = 100,      # raise to 1000 for the final run
  seed       = 1,
  outdir     = OUTDIR,
  n_cores    = 1,        # raise to parallelise bootstrap refits
  threads    = 1
)

print(result)

saveRDS(result, file.path(OUTDIR, paste0("kpower_result_", THRESHOLD, ".rds")))
write.csv(result$empirical, file.path(OUTDIR, paste0("empirical_ic_", THRESHOLD, ".csv")), row.names = FALSE)
ggplot2::ggsave(
  file.path(FIGDIR, paste0("kpower_", THRESHOLD, ".pdf")),
  result$plot, width = 7, height = 5
)
