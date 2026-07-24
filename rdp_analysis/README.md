# Recombination detection in HIV M-group subtypes with OpenRDP

This folder prepares an explicit recombination-detection analysis of the same
HIV subtype consensus alignments used for the MAST work. The MAST model detects
recombination indirectly, by showing that different sites favour different trees.
OpenRDP complements this. It scans the alignment for individual recombination
events, and for each event it reports the recombinant sequence, the two likely
parents, and the breakpoint positions. Together the two approaches answer
different questions about the same data: MAST asks how many distinct tree
histories the sites need, whereas OpenRDP asks which sequences recombined and
where.

Nothing in this folder has been run. OpenRDP is not installed in this
environment, and no recombination tool of any kind is present here apart from
MAST. The files below set up the inputs, the parameters, and the run command so
that the analysis is ready once OpenRDP is installed.

## Inputs

The analysis uses the three gap-trimmed alignments in `../data`, the same
inputs as the kpower analysis:

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

## Tool: OpenRDP

OpenRDP (https://github.com/PoonLab/OpenRDP) is a cross-platform Python
reimplementation of the method suite made popular by RDP. It runs from the
command line, which suits the existing IQ-TREE and R workflow and keeps the
analysis reproducible. It implements seven methods: RDP, GENECONV, Bootscan,
MaxChi, Chimaera, SiScan, and 3Seq. Reporting an event under several independent
methods is the usual standard of evidence in this literature, so we run all
seven and look for agreement.

### Install

OpenRDP is not present here. Install it before running, for example:

```
python3 -m pip install openrdp
```

Alternatively, clone the repository and install from source:

```
git clone https://github.com/PoonLab/OpenRDP
cd OpenRDP && python3 -m pip install --user .
```

### Run

From this folder:

```
cd mast_hiv/rdp_analysis
zsh scripts/run_openrdp.sh
```

The script scans each of the three alignments with all seven methods and writes
one comma-separated results file per threshold to `results/<threshold>/`. It
checks that OpenRDP is on the path first and reports a clear error if it is not.
The equivalent single call for one alignment is:

```
openrdp ../data/0.1_trimmed_alns_HIV_2021_subset.fasta \
  -c config/hiv_rdp.ini -o results/0.1/openrdp_0.1.csv \
  -m rdp geneconv bootscan maxchi chimaera siscan threeseq
```

## Configuration

`config/hiv_rdp.ini` holds the per-method settings. The values start from the
OpenRDP defaults, which are a reasonable choice for these moderately divergent
alignments. Two settings deserve attention when reading the results. First, the
window sizes control the resolution of breakpoint detection, and the shorter
0.9 alignment may need a smaller `win_size` if a window spans too few variable
sites. Second, the 0.9 alignment still contains gap columns, and the
`strip_gaps` and `indels_as_polymorphisms` settings decide how each method
treats them. Keep these settings identical across the three runs so that any
difference in the results reflects the data rather than the parameters.

## Outputs and interpretation

Each run writes a table with one row per detected event and the columns Method,
Start, End, Recombinant, Parent1, Parent2, and Pvalue. Start and End are the
inferred breakpoint positions in alignment coordinates. An event supported by
several methods, at consistent breakpoints, is far stronger evidence than an
event flagged by a single method. Because the HXB2 reference is the first
sequence and carries standard genome coordinates, breakpoints can be mapped back
to named HIV genome regions for interpretation.

## Relationship to the MAST and kpower analyses

The MAST analysis in `../output` found that six trees best describe the full
alignment, which points to several distinct evolutionary histories across the
genome. The kpower analysis in `../kpower_analysis` tests whether each trimmed
alignment carries enough signal to recover that number of trees. OpenRDP adds
the third piece: the specific recombinant sequences and breakpoints that give
rise to those multiple histories. Running OpenRDP on the same three alignments
lets us check whether the events it finds line up with the tree-mixture signal
that MAST detects, and whether both signals survive the same gap trimming.
