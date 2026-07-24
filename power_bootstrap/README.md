# Parametric bootstrap power test for MAST tree-number selection

This analysis measures how reliably the core MAST analysis recovers its own
selected number of trees. The core analysis fit MAST models with k = 2 to 35
trees to each gap-trimmed alignment and chose the number of trees that minimised
the Bayesian information criterion (BIC). That number differs by threshold: four
trees for 0.1, six for 0.5, and five for 0.9. The question here is whether that
choice is trustworthy. If the data really followed the selected model, would the
same procedure pick the same number of trees again? The parametric bootstrap
answers this directly by simulating data under the selected model and re-running
the selection.

The test is built to match the core analysis exactly. It reuses the core fitted
models, the core candidate trees, and the identical IQ-TREE call. Only the input
alignment changes, from the real data to a simulated replicate.

## What the test does

For each threshold t with core-selected best K\*:

1. **Read the fitted K\* model** from `output/<t>/mast_model_K*/`. This supplies
   the linked GTR substitution rates and base frequencies, the FreeRate rate
   categories (R4), the K\* tree weights, and the K\* trees with their fitted
   branch lengths.
2. **Simulate B replicate alignments** under that model with IQ-TREE's AliSim.
   Each site is assigned to one of the K\* tree classes in proportion to the
   fitted weights, then evolved along that class tree under the shared GTR+FO+R4
   process. Each replicate has the same length as the real alignment.
3. **Refit MAST for k = 2 to 10** on every replicate. The refit reuses the core
   candidate-tree subsets in `generated/<t>/tree_subsets/modeltrees_final_k.newick`
   and the identical call `iqtree3 -m "GTR+FO+R4+T" -te <subset> -wspmr -wslmr`.
   The topologies are fixed and the branch lengths are re-optimised, as in the
   core sweep.
4. **Record the BIC-minimising k** for each replicate. The power to recover K\*
   is the fraction of replicates whose minimum-BIC k equals K\*.

The refit range stops at k = 10 rather than 35. The BIC-vs-k curve in every core
sweep has a single interior minimum well below 10, so this range captures the
selected k with margin at a fraction of the cost.

## Why this replaces the earlier kpower attempt

An earlier attempt used the kpower package directly. That approach rebuilt its
own candidate trees from alignment windows and selected its own rate model, so
its K did not correspond to the core analysis. This test removes that gap. It
takes the core trees, the core model, and the core call as fixed inputs, so its
K is the same quantity the core analysis reports.

## Requirements

- IQ-TREE 3 at the path set in `power_bootstrap.R` (`PB_IQTREE` overrides it).
- R packages `here`, `ape`, and `parallel`.
- The core outputs already present in the repository: the best-K model files
  under `output/<t>/mast_model_K*/` and the candidate trees under
  `generated/<t>/tree_subsets/`.

## Running

From the repository root, so that `here::i_am` anchors on `mast_hiv.Rproj`:

```
zsh power_bootstrap/run_power_bootstrap.sh
```

This runs all three thresholds and writes per-threshold logs to
`results/logs/`. To run one threshold directly:

```
Rscript power_bootstrap/power_bootstrap.R 0.1
```

For a quick end-to-end check at reduced scale, override the defaults through
environment variables:

```
PB_B=2 PB_KMAX=4 PB_NCORES=4 Rscript power_bootstrap/power_bootstrap.R 0.1
```

Available overrides: `PB_B` (replicates), `PB_KMIN`, `PB_KMAX` (refit range),
`PB_NCORES` (concurrent fits), `PB_NT` (threads per fit), `PB_SEED`, `PB_IQTREE`.

## Running on a bigger machine

The full B = 100 analysis is heavy. On a 10-core laptop the three thresholds
take roughly three days, because each replicate needs a k = 2 to 10 sweep of
branch-length re-optimising MAST fits. A server with more cores finishes it in
proportionally less time, since the fits are independent and run in parallel.

To move it to another machine, copy the repository there, make sure IQ-TREE 3,
R, and the packages `here` and `ape` are available, then run the driver from the
repository root:

```
zsh power_bootstrap/run_power_bootstrap.sh
```

By default the script uses every detected core, one thread per fit, which gives
the best throughput for many small independent fits. Point it at the machine's
IQ-TREE and, if needed, cap the cores:

```
PB_IQTREE=/path/to/iqtree3 PB_NCORES=32 zsh power_bootstrap/run_power_bootstrap.sh
```

The pipeline is otherwise self-contained. It anchors paths with `here::i_am`, so
it must run from the repository root, and it resolves IQ-TREE in this order: the
`PB_IQTREE` variable, then `iqtree3` or `iqtree2` on `PATH`, then the local macOS
build used during development. The only inputs it reads are the core outputs
already in the repository (`output/<t>/mast_model_K*/` and
`generated/<t>/tree_subsets/`), so a fresh copy plus IQ-TREE 3 and the two R
packages is enough to run it.

## Outputs

Per threshold, written to `results/<t>/`:

| File | Contents |
|------|----------|
| `bic_by_rep_<t>.csv` | BIC for every replicate and every k |
| `recovered_K_<t>.csv` | BIC-minimising k per replicate |
| `power_summary_<t>.txt` | power to recover K\* and the distribution of recovered K |
| `recovered_K_<t>.pdf` | bar plot of recovered K across replicates |

## Interpretation

A high power, near 100 percent, means the core analysis reliably identifies its
number of trees, so the selected K is well supported. A lower power means BIC
often selects a different number of trees on data generated by the selected
model, which is a caution against reading K\* too strictly. Comparing the three
thresholds shows whether trimming gappy columns strengthens or weakens the
signal for the number of distinct tree histories.
