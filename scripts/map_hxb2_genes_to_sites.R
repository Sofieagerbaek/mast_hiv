library(ape)
library(seqinr)

# Reproduce trimCols() exactly as in scripts/make_ML_trees.Rmd
trimCols <- function(al, prop, codon = TRUE) {
  mat <- as.character(as.matrix(al))
  ntax <- nrow(mat)
  propthres <- 1 - prop
  compliantSites <- apply(mat, 2, function(x) {
    x <- as.character(x)
    compl <- (length(which(x %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) / ntax) < propthres
    return(compl)
  })
  if (codon) {
    codIDs <- rep(1:(length(compliantSites) / 3), each = 3)
    codsToKeep <- rep(FALSE, length(compliantSites))
    for (i in 1:max(codIDs)) {
      if (all(compliantSites[which(codIDs == i)])) codsToKeep[which(codIDs == i)] <- TRUE
    }
    list(keep = as.logical(codsToKeep))
  } else {
    list(keep = as.logical(compliantSites))
  }
}

genome_alignment <- read.alignment("data/HIV1_CON_2021_genome_DNA_subset.fasta", format = "fasta")
alignment_matrix <- as.DNAbin(genome_alignment)
mat <- as.character(as.matrix(alignment_matrix))
cat("Untrimmed alignment: ", nrow(mat), "taxa x", ncol(mat), "columns\n")

hxb2_row <- mat["B.FR.83.HXBLAI_IIIB_BRU.K03455", ]
is_gap <- hxb2_row %in% c("-", "n", "N", "?")
hxb2_ungapped_len <- sum(!is_gap)
cat("HXB2 ungapped length in this alignment:", hxb2_ungapped_len, "bp (should be 9719)\n")

# Map: HXB2 raw position (1..N) -> alignment column index (1..ncol, untrimmed)
col_of_hxb2_pos <- which(!is_gap)

# Gene coordinates in raw HXB2 numbering (K03455), from NCBI GenBank CDS features.
# vpu is not annotated as a CDS in the historical K03455 record (defective start codon
# in this isolate), so its range is taken from the widely-cited HXB2 landmark numbering
# used throughout the HIV literature (e.g. Korber et al. 1998, "Numbering Positions in
# HIV Relative to HXB2CG").
genes <- list(
  gag = c(790, 2292),
  pol = c(2085, 5096),   # 2085 = conventional frameshift-junction start; GenBank's own
                          # CDS begins at 2358 (post-frameshift ORF) - see pol_genbank_cds row
  vif = c(5041, 5619),
  vpr = c(5559, 5850),   # HXB2's own CDS is frameshifted/truncated (5559-5795); 5850 is
                          # the standard full-length landmark used here
  tat = c(5831, 8469),   # two exons (5831-6045, 8379-8469); the intron (6046-8378) is
                          # not part of the mature transcript
  rev = c(5970, 8653),   # two exons (5970-6045, 8379-8653); same intron caveat as tat
  vpu = c(6062, 6310),
  env = c(6225, 8795),
  nef = c(8797, 9417)    # HXB2's own copy has a premature stop at 9168; 9417 is the
                          # standard full-length landmark used here
)
pol_genbank <- c(2358, 5096)

# tat/rev individual exons (the coding parts only, excluding the intron between them)
exons <- list(
  tat_exon1 = c(5831, 6045),
  tat_exon2 = c(8379, 8469),
  rev_exon1 = c(5970, 6045),
  rev_exon2 = c(8379, 8653)
)

hxb2_max_pos <- length(col_of_hxb2_pos)

get_col_range <- function(pos_range) {
  s <- min(pos_range, hxb2_max_pos)
  e <- min(max(pos_range), hxb2_max_pos)
  c(col_of_hxb2_pos[s], col_of_hxb2_pos[e])
}

# Per-threshold: which original columns survive trimCols(), and their new (trimmed) site index
thresholds <- c("0.1", "0.5", "0.9")
keep_list <- list()
site_index_list <- list()
for (t in thresholds) {
  keep <- trimCols(alignment_matrix, as.numeric(t), codon = FALSE)$keep
  keep_list[[t]] <- keep
  site_idx <- rep(NA_integer_, length(keep))
  site_idx[keep] <- seq_len(sum(keep))
  site_index_list[[t]] <- site_idx
  cat(sprintf("threshold %s: retained %d / %d columns\n", t, sum(keep), length(keep)))
}

notes <- c(
  gag = "",
  pol = "2085 = conventional frameshift-junction start (commonly cited pol start); GenBank's own CDS begins at 2358 (post-frameshift ORF)",
  vif = "annotated as \"sor\" in the original GenBank record",
  vpr = "HXB2's own CDS is frameshifted/truncated (5559-5795); 5850 is the standard full-length landmark used here",
  tat = "two exons (5831-6045, 8379-8469), spliced; the intron (6046-8378) is not part of the mature transcript, so the mRNA is shorter than this span",
  rev = "two exons (5970-6045, 8379-8653), spliced; same intron caveat as tat",
  vpu = "not annotated as a CDS in the historical K03455 GenBank record (defective start codon in this isolate); coordinates from standard HIV literature (e.g. Korber et al. 1998)",
  env = "",
  nef = "HXB2's own copy has a premature stop at 9168; 9417 is the standard full-length landmark used here"
)

rows <- lapply(names(genes), function(g) {
  colr <- get_col_range(genes[[g]])
  row <- list(
    gene = g,
    hxb2_start = genes[[g]][1],
    hxb2_end = genes[[g]][2],
    aligned_col_start = colr[1],
    aligned_col_end = colr[2]
  )
  for (t in thresholds) {
    keep <- keep_list[[t]]; site_idx <- site_index_list[[t]]
    cols_in_span <- colr[1]:colr[2]
    kept_in_span <- cols_in_span[keep[cols_in_span]]
    row[[paste0("site_", t, "_start")]] <- site_idx[min(kept_in_span)]
    row[[paste0("site_", t, "_end")]]   <- site_idx[max(kept_in_span)]
    row[[paste0("aligned_cols_trimmed_", t)]] <- length(cols_in_span) - length(kept_in_span)
  }
  row$note <- notes[g]
  as.data.frame(row, stringsAsFactors = FALSE)
})
gene_table <- do.call(rbind, rows)

pol_row <- data.frame(
  gene = "pol_genbank_cds",
  hxb2_start = pol_genbank[1], hxb2_end = pol_genbank[2],
  aligned_col_start = get_col_range(pol_genbank)[1], aligned_col_end = get_col_range(pol_genbank)[2],
  site_0.1_start = NA, site_0.1_end = NA, aligned_cols_trimmed_0.1 = NA,
  site_0.5_start = NA, site_0.5_end = NA, aligned_cols_trimmed_0.5 = NA,
  site_0.9_start = NA, site_0.9_end = NA, aligned_cols_trimmed_0.9 = NA,
  note = "Alternative pol range using GenBank's own post-frameshift CDS start (2358) instead of the conventional frameshift-junction landmark (2085)"
)
for (t in thresholds) {
  keep <- keep_list[[t]]; site_idx <- site_index_list[[t]]
  colr <- get_col_range(pol_genbank)
  cols_in_span <- colr[1]:colr[2]
  kept_in_span <- cols_in_span[keep[cols_in_span]]
  pol_row[[paste0("site_", t, "_start")]] <- site_idx[min(kept_in_span)]
  pol_row[[paste0("site_", t, "_end")]]   <- site_idx[max(kept_in_span)]
  pol_row[[paste0("aligned_cols_trimmed_", t)]] <- length(cols_in_span) - length(kept_in_span)
}
gene_table <- rbind(gene_table, pol_row)

# tat/rev exon-only rows (excludes the intron, unlike the combined tat/rev rows above)
exon_notes <- c(
  tat_exon1 = "tat exon 1 (coding only)",
  tat_exon2 = "tat exon 2 (coding only)",
  rev_exon1 = "rev exon 1 (coding only)",
  rev_exon2 = "rev exon 2 (coding only)"
)
exon_rows <- lapply(names(exons), function(g) {
  colr <- get_col_range(exons[[g]])
  row <- list(
    gene = g,
    hxb2_start = exons[[g]][1],
    hxb2_end = exons[[g]][2],
    aligned_col_start = colr[1],
    aligned_col_end = colr[2]
  )
  for (t in thresholds) {
    keep <- keep_list[[t]]; site_idx <- site_index_list[[t]]
    cols_in_span <- colr[1]:colr[2]
    kept_in_span <- cols_in_span[keep[cols_in_span]]
    row[[paste0("site_", t, "_start")]] <- site_idx[min(kept_in_span)]
    row[[paste0("site_", t, "_end")]]   <- site_idx[max(kept_in_span)]
    row[[paste0("aligned_cols_trimmed_", t)]] <- length(cols_in_span) - length(kept_in_span)
  }
  row$note <- exon_notes[g]
  as.data.frame(row, stringsAsFactors = FALSE)
})
gene_table <- rbind(gene_table, do.call(rbind, exon_rows))

out_path <- "data/hxb2_gene_positions.csv"
write.csv(gene_table, out_path, row.names = FALSE)
cat("\nWrote:", out_path, "\n")
