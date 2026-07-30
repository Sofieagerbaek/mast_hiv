# Recombination detection in HIV M-group subtypes with OpenRDP

This folder holds an explicit recombination-detection analysis of the same HIV
subtype consensus alignments used for the MAST work. The MAST model detects
recombination indirectly, by showing that different sites favour different trees.
OpenRDP complements this. It scans the alignment for individual recombination
events, and for each event it reports the recombinant sequence, the two likely
parents, and the breakpoint positions. Together the two approaches answer
different questions about the same data: MAST asks how many distinct tree
histories the sites need, whereas OpenRDP asks which sequences recombined and
where.

**Status: complete.** All seven methods have been run on all three thresholds.
Read the "Results" section before using any of it — most of the methods do not
discriminate on this dataset, and one of them reports invalid p-values.

## Results

Scored against the alignment's own labels, treating the 20 recombinant forms
(AE, AG, 01B, ...) as true positives and the 17 pure subtypes plus the CPZ
outgroup as true negatives. Youden's J = sensitivity + specificity - 1, so
J = 0 is no better than flagging everything. Full table in
`results/method_scorecard.csv`.

| Method | Events (0.1 / 0.5 / 0.9) | Youden's J | Verdict |
|---|---|---|---|
| **3Seq** | 34 / 33 / 33 | **0.12 / 0.17 / 0.28** | **the only usable method here** |
| SiScan | 1 / 4 / 0 | 0.05 / 0.04 / n/a | perfectly specific but near-zero sensitivity; degenerate spans |
| Bootscan | 47,564 / 35,702 / 24,622 | 0 / 0 / 0 | flags all 37 sequences, no discrimination |
| RDP | 67,636 / 63,442 / 57,958 | 0 / 0 / 0 | flags all 37; **invalid p-values, see below** |
| MaxChi | 4,240 / 3,151 / 75 | 0 / 0 / -0.54 | flags all 37 at 0.1 and 0.5 |
| GENECONV | 80 / 81 / 80 | -0.21 | selective but worse than chance |
| Chimaera | 0 / 0 / 0 | n/a | zero events on all three, even after patching |

**Use 3Seq.** It is the only method with positive discrimination, and it is best
on the most aggressively trimmed alignment (0.9). Treat everything else as
unusable on this dataset rather than as weak corroboration: a method that flags
all 37 sequences, including the pure HXB2 reference, agrees with 3Seq by
construction and adds no evidence.

Two further cautions:

- **The RDP method's `Pvalue` column is not a probability.** 65,972 of its 67,636
  events at threshold 0.1 fall outside [0, 1] — values like 4743029.99. This is
  an upstream OpenRDP defect, not a setup problem: `openrdp/rdp.py` carries its
  own `# TODO check for valid p-val ... output has given values greater than 31
  before`, and the `max_pvalue` setting under `[RDP]` in `config/hiv_rdp.ini` is
  documented in the source as unused. The `pval_out_of_range` column in the
  scorecard counts these per method; RDP is the only one affected.
- **Breakpoint positions are more trustworthy than sequence assignments.**
  OpenRDP's choice of which sequence in a triplet is the recombinant and which
  are the parents is unreliable here. Prefer the `Start`/`End` columns over
  `Recombinant`/`Parent1`/`Parent2`. SiScan is an extreme case: it returns
  `End = Start - 1` (e.g. 9540 → 9539), so it yields no usable interval at all.

The OpenRDP README itself warns the tool is "probably not yet ready for
practical use." That is consistent with what we see.

## Inputs

The analysis uses the three gap-trimmed alignments in `../data`, the same inputs
as the MAST core analysis and the power bootstrap in `../power_bootstrap`:

| Threshold | File | Aligned length |
|-----------|------|----------------|
| 0.1 | `../data/0.1_trimmed_alns_HIV_2021_subset.fasta` | 9,546 sites |
| 0.5 | `../data/0.5_trimmed_alns_HIV_2021_subset.fasta` | 8,941 sites |
| 0.9 | `../data/0.9_trimmed_alns_HIV_2021_subset.fasta` | 8,546 sites |

Each alignment holds the same 37 sequences: the HXB2 reference
(`B.FR.83.HXBLAI_IIIB_BRU.K03455`) followed by 36 subtype and circulating
recombinant form consensus sequences. The sequence names are short and free of
spaces and special characters, so OpenRDP reads them without renaming. The three
files differ only in how aggressively gappy columns were removed, which lets us
ask whether the recombination signal is robust to alignment trimming.

## Install

OpenRDP (https://github.com/PoonLab/OpenRDP) is a cross-platform Python
reimplementation of the method suite made popular by RDP. It implements the seven
methods above.

**It is not on PyPI** — `pip install openrdp` fails. Install from source:

```
git clone https://github.com/PoonLab/OpenRDP
cd OpenRDP && python3 -m pip install --user .
```

Its two compiled dependencies, 3Seq and GENECONV, ship as prebuilt native Linux
binaries inside the package, so nothing needs compiling.

### Required local patches

**The upstream release has four bugs that materially affect these results.** They
are fixed by `patches/openrdp-local-fixes.patch` (against openrdp 0.1.0), which
must be applied to the installed package before running anything:

The patch paths are `a/openrdp/<file>.py`, so apply it with `-p1` from the
directory *containing* the installed `openrdp` package, not from inside it:

```
python3 -m pip install --user .            # from the OpenRDP clone
cd "$(python3 -c 'import openrdp, os; print(os.path.dirname(os.path.dirname(openrdp.__file__)))')"
patch -p1 < /path/to/mast_hiv/rdp_analysis/patches/openrdp-local-fixes.patch
```

Each patched site is marked `# LOCAL PATCH`. Note that `pip install --upgrade
openrdp` silently discards all of it. What the patch fixes:

1. **`common.calculate_chi2`** — the validity guard tested the wrong condition.
   It let through tables that make `scipy.chi2_contingency` raise, which crashed
   MaxChi on all three thresholds, while *discarding* perfectly-associated tables
   like `[[5,0],[0,3]]` — the strongest MaxChi/Chimaera signal there is. Now
   tests row and column sums, which is the actual requirement.
2. **`maxchi.py`** and 3. **`chimaera.py`**, at the `left_peak`/`right_peak`
   lookup after greedy window expansion — the expand guard tests
   `right < len(...)` and *then* increments, so `right` can index one past the
   end while `left` can fall below 0 and wrap around from the back. The result
   was either an `IndexError` (which is what crashed MaxChi at 0.5 and 0.9) or a
   silently wrong value. Both indices are now clamped to the valid range, which
   also keeps the derived breakpoint coordinates inside the alignment.
4. **`siscan.find_signal`** — an infinite loop that leaked roughly 5 GB/h and got
   the process OOM-killed. `find_interval` returns `end - 1`, and when it breaks
   on its first iteration (whenever `close[ind][ps] == 0`, since 0 is falsy) it
   returns `ind - 1`; the old `ind = new` then walked backwards, the following
   `ind += 1` restored the original index, and the enclosing `while` re-tested
   the same still-true condition, appending to `m1_intervals` forever. Now
   `ind = max(new, ind)`, with the same fix on the `ind2`/maj2 loop.

## Run

Each `openrdp` invocation is single-core. Its `np` setting (under `[Bootscan]`)
is read but never wired to any parallel code path, and its only real parallelism
hook is `mpi4py`, which falls back to serial without warning when absent. So
rather than one core working through 7 methods per threshold, the driver launches
all 3 thresholds × 7 methods = 21 independent single-core jobs at once:

```
cd rdp_analysis
bash run_openrdp_parallel.sh     # 21 jobs -> results/<t>/parts/, then merges
```

It writes one CSV per method to `results/<t>/parts/` and then merges each
threshold's parts into a combined `results/<t>/openrdp_<t>.csv`. Expect this to
take the better part of a day; a single-method smoke test (RDP on the smallest
0.9 alignment) takes ~14 minutes on one core, and Bootscan
(`num_replicates=100`) and SiScan (`pvalue_perm_num=1000`) are far slower.

To rebuild the combined CSVs from whatever parts exist, at any time:

```
bash merge_all.sh
```

It reports which methods contributed and which are missing, so a partial
analysis stays visible rather than looking complete.

### SiScan needs chunking

Even with the infinite-loop fix, SiScan on 7,770 triplets is too slow to finish
in one process. `rerun_siscan.sh` splits a threshold across 9 concurrent
`siscan_chunk.py` processes and `merge_siscan_chunks.py` reassembles them,
refusing to write a partial merge. Chunking is exact rather than approximate:
`Siscan.execute` reseeds per triplet, so a triplet's result never depends on how
many triplets ran before it. Each chunk took ~8.4 h.

`rerun_failed.sh` re-runs just the methods that produced no usable output in a
first pass. Both rerun scripts capture exit statuses correctly, which the
original driver did not — it wrote `echo "... $(date) ... (exit $?)"`, and the
command substitution resets `$?` before it is expanded, so every job was logged
as "exit 0" including the ones that crashed or were OOM-killed. **Check part-file
row counts, not the driver log's exit codes.**

## Figures

```
Rscript rdp_analysis/summarise_and_plot.R    # from the repo root
```

Writes `results/method_scorecard.csv` and five PDFs to `results/figures/`:
breakpoint density by method, a breakpoint histogram for the selective methods
only, a mosaic map of events along the genome, a method-agreement plot, and the
scorecard. It reads `results/<t>/parts/`, one file per method, not the combined
CSVs — see the note below on why both are committed. It also reads the
untrimmed `../data/HIV1_CON_2021_genome_DNA_subset.fasta` to map trimmed
coordinates back to genome positions.

Requires R with `here`, `readr`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, and
`scales`.

## Configuration

`config/hiv_rdp.ini` holds the per-method settings, starting from the OpenRDP
defaults, which are reasonable for these moderately divergent alignments. Two
settings deserve attention when reading the results. First, the window sizes
control the resolution of breakpoint detection, and the shorter 0.9 alignment may
need a smaller `win_size` if a window spans too few variable sites. Second, the
0.9 alignment still contains gap columns, and `strip_gaps` and
`indels_as_polymorphisms` decide how each method treats them. Keep these settings
identical across the three runs so that any difference in the results reflects
the data rather than the parameters. As noted above, `max_pvalue` under `[RDP]`
is inert.

## What is in `results/`

| Path | Contents |
|---|---|
| `<t>/openrdp_<t>.csv` | all methods merged, one row per event: `Method,Start,End,Recombinant,Parent1,Parent2,Pvalue` |
| `<t>/parts/openrdp_<t>_<method>.csv` | the same rows split per method, as each method emitted them |
| `method_scorecard.csv` | per method and threshold: event count, out-of-range p-values, median span, sensitivity, specificity, Youden's J |
| `figures/*.pdf` | the five plots from `summarise_and_plot.R` |
| `logs/parallel_driver.log`, `logs/siscan_rerun2/driver.log` | which jobs ran and when |

`parts/` is not redundant with the combined CSVs, despite holding the same rows.
A header-only part file means "this method ran and found nothing" (Chimaera on
all three thresholds), whereas an absent part file means "this method never
produced output" (MaxChi before patching). The combined CSV cannot express that
difference, and `summarise_and_plot.R` relies on it. The per-job run logs are
gitignored — 23 MB of mostly RDP and Bootscan progress spam — but the two driver
logs are kept.

## Relationship to the MAST and power-bootstrap analyses

The MAST analysis in `../output` selects 4, 6, and 5 trees for thresholds 0.1,
0.5, and 0.9 respectively, which points to several distinct evolutionary
histories across the genome. The parametric bootstrap in `../power_bootstrap`
confirms that those counts are recoverable (98%, 100%, 100%). OpenRDP was
intended to supply the third piece: the specific recombinant sequences and
breakpoints giving rise to those histories.

That third piece is weaker than hoped. Only 3Seq discriminates at all, so the
event-level picture rests on a single method rather than on the cross-method
agreement that this literature normally expects. Treat the OpenRDP results as
loose corroboration that recombination is present and roughly where, not as an
independent confirmation of the MAST tree count.
