# Examining recombination in HIV M-group subtypes with IQ-tree's mixture across sites and trees (MAST) model

Scripts and files used for the analysis of a subset of 37 of the 2021 HIV subtype 1 consensus sequences (https://www.hiv.lanl.gov/content/index) using mixture models.

The repository holds three analyses:

| Folder | Analysis | Status |
|---|---|---|
| `scripts/`, `output/`, `generated/` | the core MAST analysis: how many distinct tree histories the sites need | complete |
| `power_bootstrap/` | parametric bootstrap: is that number of trees reliably recoverable? | complete |
| `rdp_analysis/` | OpenRDP: which sequences recombined, and where? | complete, with caveats |

Each of the two newer folders has its own README with results, requirements, and
how to re-run it. Start there rather than here for those.

#### Layout

- Scripts for the core analysis in the `scripts/` folder
- Output from the full MAST model run on all 35 initial trees, and MAST output for every k from 2 to 35, in the `output/` folder
- Intermediate trees and candidate tree subsets in the `generated/` folder
- Alignments in `data/`

#### Core analysis steps

Trimming of the alignment and tree building was done with the `make_ML_trees.Rmd` script. Tree-building output for the 18 initial and 17 shifted trees is under `generated/<threshold>/18_ML_trees/` and `17_ML_trees/`, and the 35 combined trees are in `generated/<threshold>/35_trimmed_ML_trees.newick`.

IQ-tree's MAST (http://www.iqtree.org/doc/Complex-Models) model was run on those 35 newick trees and the trimmed alignment in the command line with the call:

    iqtree2 -s $FASTA -m "GTR+FO+R4+T" -te $TREES --prefix $OUT1 -fast -nt AUTO -wspmr -wslmr -redo

The MAST analysis was rerun on sets of trees from k=2 to k=34 iteratively adding trees in order of their weight. Bayesian Information Criterion (BIC) from these was used to determine the optimal number of trees.

The analysis was run separately on three alignments that differ only in how aggressively gappy columns were trimmed, and the selected number of trees differs between them:

| Threshold | Alignment | Aligned length | Selected K\* | MAST output |
|---|---|---|---|---|
| 0.1 | `data/0.1_trimmed_alns_HIV_2021_subset.fasta` | 9,546 sites | 4 | `output/0.1/mast_model_4/` |
| 0.5 | `data/0.5_trimmed_alns_HIV_2021_subset.fasta` | 8,941 sites | 6 | `output/0.5_output/mast_model_6/` |
| 0.9 | `data/0.9_trimmed_alns_HIV_2021_subset.fasta` | 8,546 sites | 5 | `output/0.9/mast_model_5/` |

Note the inconsistent `output/0.5_output/` folder name, which breaks the pattern of the other two. Scripts that read these paths handle it with an explicit lookup rather than by string interpolation (see `power_bootstrap/power_bootstrap.R`), so keep that in mind when adding a threshold.

Per-site results for each selected model — the posterior tree-class assignment per site and the per-site log-likelihoods — are in the `.siteprob` and `.sitelh` files in those folders, and are the basis for the across-sites figures.

Post-processing and most visualisations are performed in `scripts/newer_MAST_post_processing.Rmd`; the knitted HTML is also provided. An earlier version is kept as `scripts/2503_MAST_post_processing.Rmd`.

Some other smaller visualisation scripts are present in the scripts folder, such as information criteria and weight plotting — as well as the shell script used to run iqtree MAST runs (first and second) and scf calculations.

#### Follow-up analyses

**`power_bootstrap/`** — tests whether the core analysis recovers its own selected K\*, by simulating under the selected model and re-running the BIC selection. Power is 98% at threshold 0.1 and 100% at 0.5 and 0.9, so the selected tree counts are well supported. See `power_bootstrap/README.md`.

**`rdp_analysis/`** — runs OpenRDP's seven recombination-detection methods on the same three alignments to identify individual events and breakpoints. **Read `rdp_analysis/README.md` before using these results.** Only 3Seq discriminates between the known recombinant forms and the pure subtypes on this dataset; Bootscan, RDP and MaxChi flag all 37 sequences indiscriminately, and the RDP method's reported p-values are not valid probabilities. The analysis also requires four local bug-fix patches to OpenRDP, supplied in `rdp_analysis/patches/`.
