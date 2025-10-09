here::i_am("scripts/ML_tree_compar.R")
library(here)
library(dplyr)
library(phangorn)
library(ggplot2)
library(reshape2)
library(viridis)

# Load trees
ml18 <- read.tree("data/18_ML_trees/all_18_trimmed_ML_trees.newick")
ml17 <- read.tree("data/17_ML_trees/17_shifted_windows_trimmed_ML_trees.newick")

# Calculate distances
rf_distm <- matrix(NA, 18, 17)
for(i in 1:18) for(j in 1:17) rf_distm[i, j] <- RF.dist(ml18[[i]], ml17[[j]], normalize = T)

# Plot
rf_dist <- melt(rf_distm)
names(rf_dist) <- c("original", "shifted", "RFscore_norm")
rf_dist <- rf_dist %>%
  mutate(abs_shift_dist = abs(shifted - original),
         shift_dist_d = as.factor(ifelse(abs_shift_dist == 0, "0",
                                      ifelse(abs_shift_dist == 1, paste0("\U00B1","1"),
                                             paste0("\U003E","\U7C","\U00B1","1","\U7C")))),
         shift_dist_d = factor(shift_dist_d,
                               levels=c("0", paste0("\U00B1","1"), paste0("\U003E","\U7C","\U00B1","1","\U7C"))))
# Heatmap
rf_heatmap <- ggplot(rf_dist, aes(y=original, x=shifted, fill=RFscore_norm)) + 
  geom_tile() + 
  ggtitle("Normalised RF score\nML 18 vs ML shifted 17 ") +
  scale_fill_viridis(option = 3) +
  theme_classic() +
  theme(legend.position = "bottom") +
  coord_fixed()

# Histogram
rf_density <- ggplot(data = rf_dist, mapping = aes(x = RFscore_norm, fill = shift_dist_d)) + 
  geom_density(position = "identity", alpha = 0.5) +
  ggtitle("Density by shift category") +
  scale_fill_viridis_d(option = 7) +
  theme_classic() +
  theme(legend.position = "bottom") 

# Scatterplot
rf_shift <- ggplot(data = rf_dist, mapping = aes(y = RFscore_norm, x = abs_shift_dist, colour = shift_dist_d)) + 
  geom_point() +
  ggtitle("Normalised RF score\nvs shift distance") +
  scale_colour_viridis_d(option = 7) +
  theme_classic() +
  theme(legend.position = "bottom") 

p <- patchwork::wrap_plots(list(rf_heatmap, rf_shift, rf_density))

# Save PDF coz the ggsave can't handle character spacing properly
savePDF <- function (p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(p)
  invisible(grDevices::dev.off())
}

savePDF(p, here("figures/ml_tree_compar.pdf"), width = 10, height = 5)
