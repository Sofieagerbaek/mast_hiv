#!/usr/bin/env Rscript
# Summarise the OpenRDP output and draw breakpoint maps along the HIV-1 genome.
#
# Two jobs:
#
#  1. Score each method against the labels the alignment already carries. The 37
#     sequences are subtype consensus sequences, so their names say whether they
#     are a pure subtype (A, B, C, ... ) or a recombinant form (AE = CRF01_AE,
#     AG = CRF02_AG, 01B, A2D, ...). That gives a usable, if imperfect, truth set
#     for asking which methods discriminate and which just flag everything.
#
#  2. Plot where breakpoints fall along the genome. Positions reported by OpenRDP
#     are columns of the trimmed alignment, which differ between the three
#     thresholds and mean nothing biologically, so everything is converted to
#     HXB2 nucleotide coordinates first (the HXB2 reference is the first sequence
#     in each alignment). That makes the three thresholds comparable to each
#     other and lets the plots carry a gene track.
#
# Reads the per-method files in results/<threshold>/parts/ rather than the merged
# per-threshold CSV, so it works on a partial analysis and never races the merge.
#
# Usage: Rscript summarise_and_plot.R

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

THRESHOLDS <- c("0.1", "0.5", "0.9")
METHODS <- c("rdp", "geneconv", "bootscan", "maxchi", "chimaera", "siscan", "threeseq")
RDP_DIR <- here("rdp_analysis")
OUT_DIR <- file.path(RDP_DIR, "results", "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Ground truth from the sequence names
# ---------------------------------------------------------------------------
# Pure subtypes and sub-subtypes: a single letter optionally followed by a digit.
# Recombinant forms: anything naming two or more parents (AE, AG, BC, 01B, A2D,
# A5U, C2U, cpx, ...). CPZ is the chimpanzee SIV outgroup, not an M-group
# recombinant, and is scored with the non-recombinants.
PURE <- c("B.FR.83.HXBLAI_IIIB_BRU.K03455", "A", "A1", "A2", "A3", "A6", "B", "C",
          "D", "F1", "F2", "G", "H", "J", "K", "L", "CPZ")
RECOMB <- c("AE", "AG", "cpx", "DF", "BC", "CD", "BF", "BG", "01B", "A2D", "01A1",
            "A5U", "AD", "BF1", "02G", "A1D", "02A", "C2U", "02B", "01C")

# ---------------------------------------------------------------------------
# HXB2 gene map (K03455 coordinates)
# ---------------------------------------------------------------------------
# Three rows, because vif/vpr/tat/rev/vpu all crowd into 5.0-6.3 kb and overlap
# each other; spreading them out keeps the labels legible.
GENES <- tribble(
  ~gene,  ~start, ~end,  ~row,
  "5'LTR",     1,   634,  1,
  "gag",     790,  2292,  1,
  "vif",    5041,  5619,  1,
  "vpu",    6062,  6310,  1,
  "nef",    8797,  9417,  1,
  "pol",    2085,  5096,  2,
  "vpr",    5559,  5850,  2,
  "env",    6225,  8795,  2,
  "3'LTR",  9086,  9719,  2,
  "tat",    5831,  6045,  3,
  "rev",    5970,  6045,  3
)
N_GENE_ROWS <- 3

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- grepl("^>", lines)
  ids <- sub("^>", "", lines[hdr])
  grp <- cumsum(hdr)
  seqs <- vapply(split(lines[!hdr], grp[!hdr]),
                 function(x) paste(x, collapse = ""), character(1))
  setNames(toupper(seqs), ids)
}

# Map trimmed alignment column -> HXB2 nucleotide position.
#
# Counting non-gap characters along the trimmed HXB2 row does NOT give HXB2
# coordinates: trimming deletes columns where HXB2 itself has a base, so the
# count runs short (by 365 nt at threshold 0.1, 1,185 nt at 0.9) and a gene track
# drawn against it would be displaced by up to a kilobase. The untrimmed
# alignment carries the complete 9,719 nt HXB2 genome, so go through it.
#
# Trimming only ever removes whole columns, so each trimmed column is identical
# to some untrimmed column and the kept columns stay in order. Matching on the
# full 37-character column string, left to right, recovers which ones were kept.
column_strings <- function(aln) {
  m <- do.call(rbind, strsplit(unname(aln), ""))
  apply(m, 2, paste, collapse = "")
}

build_hxb2_maps <- function() {
  untrimmed <- read_fasta(here("data", "HIV1_CON_2021_genome_DNA_subset.fasta"))
  ucols <- column_strings(untrimmed)
  hxb2_row <- strsplit(untrimmed[[1]], "")[[1]]
  hxb2_pos <- cumsum(hxb2_row != "-" & hxb2_row != ".")
  stopifnot(max(hxb2_pos) == 9719)  # complete HXB2 genome, as expected

  maps <- list()
  for (t in THRESHOLDS) {
    trimmed <- read_fasta(here("data",
                               sprintf("%s_trimmed_alns_HIV_2021_subset.fasta", t)))
    stopifnot(identical(names(trimmed), names(untrimmed)))
    tcols <- column_strings(trimmed)

    kept <- integer(length(tcols))
    j <- 1L
    for (i in seq_along(tcols)) {
      while (j <= length(ucols) && ucols[j] != tcols[i]) j <- j + 1L
      if (j > length(ucols))
        stop(sprintf("column %d of the %s alignment has no match in the untrimmed alignment", i, t))
      kept[i] <- j
      j <- j + 1L
    }
    stopifnot(all(diff(kept) > 0))

    maps[[t]] <- hxb2_pos[kept]
    message(sprintf("  %s: %d columns -> HXB2 %d..%d",
                    t, length(tcols), min(maps[[t]]), max(maps[[t]])))
  }
  maps
}

message("Building HXB2 coordinate maps ...")
maps <- build_hxb2_maps()

# ---------------------------------------------------------------------------
# Load events
# ---------------------------------------------------------------------------
load_events <- function() {
  rows <- list()
  for (t in THRESHOLDS) {
    for (m in METHODS) {
      f <- file.path(RDP_DIR, "results", t, "parts",
                     sprintf("openrdp_%s_%s.csv", t, m))
      if (!file.exists(f) || file.info(f)$size == 0) {
        message(sprintf("  %s/%-9s no output", t, m))
        next
      }
      d <- suppressWarnings(read_csv(f, col_types = cols(.default = col_character()),
                                     progress = FALSE))
      if (nrow(d) == 0) {
        message(sprintf("  %s/%-9s 0 events", t, m))
        next
      }
      d <- d %>%
        transmute(
          threshold = t,
          method = Method,
          start = suppressWarnings(as.numeric(Start)),
          end = suppressWarnings(as.numeric(End)),
          recombinant = Recombinant,
          parent1 = Parent1,
          parent2 = Parent2,
          pvalue = suppressWarnings(as.numeric(Pvalue))
        ) %>%
        filter(!is.na(start), !is.na(end))
      message(sprintf("  %s/%-9s %d events", t, m, nrow(d)))
      rows[[length(rows) + 1]] <- d
    }
  }
  if (!length(rows)) stop("no OpenRDP output found under results/*/parts/")
  bind_rows(rows)
}

message("Loading events ...")
ev <- load_events()

# Convert to HXB2 coordinates. Clamp to the alignment width first: a couple of
# methods can report an index one past the end.
ev <- ev %>%
  rowwise() %>%
  mutate(
    ncol_aln = length(maps[[threshold]]),
    hxb2_start = maps[[threshold]][max(1, min(round(start), ncol_aln))],
    hxb2_end   = maps[[threshold]][max(1, min(round(end),   ncol_aln))]
  ) %>%
  ungroup() %>%
  mutate(
    truth = case_when(recombinant %in% RECOMB ~ "recombinant form",
                      recombinant %in% PURE ~ "pure subtype",
                      TRUE ~ "unclassified"),
    # A p-value is only usable if it is actually a probability. The RDP method's
    # is not -- OpenRDP's rdp.py carries a TODO about exactly this and emits
    # values in the millions -- so flag rather than silently filter.
    pvalue_valid = !is.na(pvalue) & pvalue >= 0 & pvalue <= 1
  )

# ---------------------------------------------------------------------------
# Method scorecard
# ---------------------------------------------------------------------------
n_recomb <- length(RECOMB)
n_pure <- length(PURE)

scorecard <- ev %>%
  filter(!pvalue_valid | pvalue <= 0.05) %>%
  group_by(threshold, method) %>%
  summarise(
    events = n(),
    pval_out_of_range = sum(!pvalue_valid),
    median_span_nt = median(hxb2_end - hxb2_start),
    flagged_recomb = n_distinct(recombinant[truth == "recombinant form"]),
    flagged_pure = n_distinct(recombinant[truth == "pure subtype"]),
    .groups = "drop"
  ) %>%
  mutate(
    sensitivity = flagged_recomb / n_recomb,
    specificity = (n_pure - flagged_pure) / n_pure,
    # Youden's J: 0 for a method that flags everything (or nothing), 1 for perfect.
    youden_j = sensitivity + specificity - 1
  ) %>%
  arrange(method, threshold)

write_csv(scorecard, file.path(RDP_DIR, "results", "method_scorecard.csv"))
message("\n=== Method scorecard (events with a usable p<=0.05, or no usable p-value) ===")
print(as.data.frame(scorecard), row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# Plot theme and gene track
# ---------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold"),
                  plot.title = element_text(face = "bold")))

# Draws the gene map in the band [ymin, ymax], below the data. Rows are stacked
# from the bottom up so the track reads like the usual HIV genome diagram.
gene_track <- function(ymin, ymax) {
  band <- (ymax - ymin) / N_GENE_ROWS
  g <- GENES %>%
    mutate(lo = ymin + (N_GENE_ROWS - row) * band + 0.08 * band,
           hi = ymin + (N_GENE_ROWS - row + 1) * band - 0.08 * band)
  list(
    geom_rect(data = g, aes(xmin = start, xmax = end, ymin = lo, ymax = hi),
              inherit.aes = FALSE, fill = "grey88", colour = "grey55",
              linewidth = 0.2),
    geom_text(data = g, aes(x = (start + end) / 2, y = (lo + hi) / 2, label = gene),
              inherit.aes = FALSE, size = 2.1, colour = "grey20")
  )
}

# Show only non-negative breaks, so the gene track's band does not put negative
# "breakpoints per 100 nt" on the axis.
positive_breaks <- function(lims) {
  b <- pretty(c(0, max(lims)))
  b[b >= 0]
}

# Gene map as vertical bands spanning the full panel height, with the names along
# the top. Unlike gene_track() this does not depend on the y scale, so it works
# under facet scales = "free_y". Only the major ORFs are banded -- shading the
# five overlapping small genes in 5.0-6.3 kb would just be a grey smear.
# nef (8797-9417) and the 3'LTR (9086-9719) overlap by 330 nt and their labels
# collide, so band nef only and let the 3' end read as nef.
BAND_GENES <- GENES %>% filter(gene %in% c("5'LTR", "gag", "pol", "vif", "vpu",
                                           "env", "nef"))
gene_bands <- function(label = TRUE) {
  out <- list(
    geom_rect(data = BAND_GENES,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.16),
    geom_vline(data = BAND_GENES, aes(xintercept = start),
               colour = "grey75", linewidth = 0.15)
  )
  if (label)
    out <- c(out, list(
      geom_text(data = BAND_GENES, aes(x = (start + end) / 2, y = Inf, label = gene),
                inherit.aes = FALSE, vjust = 1.4, size = 2, colour = "grey30")
    ))
  out
}

# Methods worth plotting: those that produced events.
plot_methods <- sort(unique(ev$method))

# ---------------------------------------------------------------------------
# Figure 1: breakpoint density along the genome, per method
# ---------------------------------------------------------------------------
# Each event contributes two breakpoints, its start and its end.
bp <- ev %>%
  filter(!pvalue_valid | pvalue <= 0.05) %>%
  select(threshold, method, hxb2_start, hxb2_end) %>%
  pivot_longer(c(hxb2_start, hxb2_end), values_to = "pos") %>%
  select(threshold, method, pos)

p1 <- ggplot(bp, aes(x = pos, colour = threshold, fill = threshold)) +
  geom_density(alpha = 0.12, adjust = 0.4, linewidth = 0.5) +
  facet_wrap(~method, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.01)) +
                     coord_cartesian(xlim = c(1, 9719)) +
  labs(title = "Where each method places recombination breakpoints",
       subtitle = "Density of breakpoint positions in HXB2 coordinates; both ends of every event counted",
       x = "HXB2 position (nt)", y = "breakpoint density",
       colour = "gap-trim\nthreshold", fill = "gap-trim\nthreshold")
ggsave(file.path(OUT_DIR, "breakpoint_density_by_method.pdf"), p1,
       width = 9, height = 1.6 * length(plot_methods) + 1.5, limitsize = FALSE)

# ---------------------------------------------------------------------------
# Figure 2: breakpoint histogram with the gene track, selective methods only
# ---------------------------------------------------------------------------
# A method that calls every one of the 37 sequences a recombinant has told us
# nothing about which ones are, and its breakpoints are spread over the whole
# genome, so a histogram of them is near-uniform and swamps everything else.
# Split the methods on exactly that: does it saturate the truth set, or not.
saturating <- scorecard %>%
  filter(flagged_pure == n_pure, flagged_recomb == n_recomb) %>%
  distinct(threshold, method)
selective <- scorecard %>%
  anti_join(saturating, by = c("threshold", "method")) %>%
  distinct(threshold, method)

message(sprintf("\nSaturating (flag all %d recombinants AND all %d pure subtypes): %s",
                n_recomb, n_pure,
                paste(sprintf("%s@%s", saturating$method, saturating$threshold),
                      collapse = ", ")))
message(sprintf("Selective: %s",
                paste(sprintf("%s@%s", selective$method, selective$threshold),
                      collapse = ", ")))

if (nrow(selective)) {
  bp_sel <- bp %>% semi_join(selective, by = c("threshold", "method"))
  p2 <- ggplot(bp_sel, aes(x = pos)) +
    gene_bands() +
    geom_histogram(aes(fill = method), binwidth = 100, colour = NA) +
    # free_y: MaxChi's single 2.4 kb spike is 4x anything else and flattens the
    # other two methods to invisibility on a shared scale.
    facet_grid(method ~ threshold, scales = "free_y") +
    scale_x_continuous(labels = comma, expand = expansion(mult = 0.01)) +
                       coord_cartesian(xlim = c(1, 9719)) +
    scale_y_continuous(breaks = positive_breaks) +
    guides(fill = "none") +
    labs(title = "Breakpoint positions, methods that discriminate at all",
         subtitle = paste("100 nt bins, HXB2 coordinates, major genes shaded, free y scales.",
                          "\nMethods that flagged all 37 sequences are excluded; see breakpoint_density_by_method.pdf for those."),
         x = "HXB2 position (nt)", y = "breakpoints per 100 nt")
  ggsave(file.path(OUT_DIR, "breakpoint_histogram_selective.pdf"), p2,
         width = 10, height = 1.8 * n_distinct(selective$method) + 2.2,
         limitsize = FALSE)
}

# ---------------------------------------------------------------------------
# Figure 3: mosaic map -- one row per sequence, bars where events were called
# ---------------------------------------------------------------------------
# This is the figure that answers "which sequence is recombinant, and over which
# part of the genome" -- the HIV CRF mosaic map idiom.
mosaic <- ev %>%
  filter(!pvalue_valid | pvalue <= 0.05, truth != "unclassified") %>%
  semi_join(selective, by = c("threshold", "method")) %>%
  mutate(recombinant = factor(recombinant,
                              levels = c(sort(intersect(RECOMB, recombinant)),
                                         sort(intersect(PURE, recombinant)))))
if (nrow(mosaic)) {
  p3 <- ggplot(mosaic, aes(y = recombinant, colour = truth)) +
    gene_bands(label = FALSE) +
    geom_segment(aes(x = hxb2_start, xend = hxb2_end,
                     yend = recombinant), linewidth = 2.4, alpha = 0.55) +
    facet_grid(. ~ method + threshold, labeller = label_wrap_gen()) +
    scale_x_continuous(labels = comma, expand = expansion(mult = 0.01)) +
                       coord_cartesian(xlim = c(1, 9719)) +
    scale_colour_manual(values = c("recombinant form" = "#2b6cb0",
                                   "pure subtype" = "#c05621")) +
    labs(title = "Detected recombinant regions per sequence",
         subtitle = "Sequences named as recombinant forms in blue, pure subtype consensus sequences in orange (calls on those are likely false positives or misassigned parents)",
         x = "HXB2 position (nt)", y = NULL, colour = "sequence is a")
  ggsave(file.path(OUT_DIR, "mosaic_map.pdf"), p3,
         width = 4 * length(selective) * length(THRESHOLDS) + 3,
         height = 0.22 * n_distinct(mosaic$recombinant) + 3, limitsize = FALSE)
}

# ---------------------------------------------------------------------------
# Figure 4: method agreement per position
# ---------------------------------------------------------------------------
# The usual standard of evidence in this literature is an event supported by
# several independent methods. Count, for each 100 nt bin, how many methods place
# at least one breakpoint there.
# Counting agreement over ALL methods is useless here: the saturating methods put
# a breakpoint in nearly every bin, so the count sits at 3-5 genome-wide and
# carries no information. Restrict to the selective methods.
agree <- bp %>%
  semi_join(selective, by = c("threshold", "method")) %>%
  mutate(bin = floor(pos / 100) * 100) %>%
  distinct(threshold, method, bin) %>%
  count(threshold, bin, name = "n_methods")

n_sel_per_thr <- selective %>% count(threshold, name = "n_available")

p4 <- ggplot(agree, aes(x = bin, y = n_methods, fill = n_methods)) +
  geom_col(width = 95) +
  geom_text(data = n_sel_per_thr, inherit.aes = FALSE,
            aes(x = 300, y = Inf, label = sprintf("%d selective methods", n_available)),
            vjust = 1.6, hjust = 0, size = 2.8, colour = "grey35") +
  facet_wrap(~threshold, ncol = 1) +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.01)) +
                     coord_cartesian(xlim = c(1, 9719)) +
  scale_y_continuous(breaks = function(l) seq(0, ceiling(max(l)))) +
  scale_fill_viridis_c(option = "mako", direction = -1, guide = "none") +
  labs(title = "How many of the discriminating methods place a breakpoint in each 100 nt bin",
       subtitle = "Saturating methods excluded -- including them puts a breakpoint in almost every bin and the count becomes uninformative",
       x = "HXB2 position (nt)", y = "number of methods")
ggsave(file.path(OUT_DIR, "method_agreement.pdf"), p4, width = 9, height = 6)

# ---------------------------------------------------------------------------
# Figure 5: the scorecard itself
# ---------------------------------------------------------------------------
sc_long <- scorecard %>%
  select(threshold, method, sensitivity, specificity, youden_j) %>%
  pivot_longer(c(sensitivity, specificity, youden_j),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         sensitivity = "sensitivity\n(recombinant forms found)",
                         specificity = "specificity\n(pure subtypes not flagged)",
                         youden_j = "Youden's J\n(sens + spec - 1)"))

p5 <- ggplot(sc_long, aes(x = value, y = method, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_fill_brewer(palette = "Blues") +
  labs(title = "How well does each method separate recombinant forms from pure subtypes?",
       subtitle = paste0("Truth set is the sequence naming: ", n_recomb,
                         " recombinant forms vs ", n_pure,
                         " pure subtypes/outgroup. Youden's J is 0 for a method that flags everything."),
       x = NULL, y = NULL, fill = "gap-trim\nthreshold")
ggsave(file.path(OUT_DIR, "method_scorecard.pdf"), p5,
       width = 11, height = 0.42 * n_distinct(sc_long$method) + 2.4, limitsize = FALSE)

message(sprintf("\nWrote figures (PDF) to %s", OUT_DIR))
message("Wrote scorecard table to results/method_scorecard.csv")
