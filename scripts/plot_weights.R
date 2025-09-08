here::i_am("scripts/plot_weights.R")
library(here)
library(dplyr)
library(magrittr)

# Load weights from .iqtree
weights_fast <- system("grep '^Tree weights' /maps/projects/racimolab/people/fnl270/GitHub/mast_hiv/output/mast_35/mast_35_fast.iqtree | sed 's/Tree weights: //; s/,//g'", intern = TRUE)
weights_normal <- system("grep '^Tree weights' /maps/projects/racimolab/people/fnl270/GitHub/mast_hiv/output/mast_35/mast_35.iqtree | sed 's/Tree weights: //; s/,//g'", intern = TRUE)
weights_fast <- as.numeric(strsplit(weights_fast, " +")[[1]])
weights_normal <- as.numeric(strsplit(weights_normal, " +")[[1]])

weights_fast <- reshape2::melt(weights_fast, value.name = "weight") %>%
    mutate(group = "fast") %>%
    mutate(tree = 1:35)

weights_normal <- reshape2::melt(weights_normal, value.name = "weight") %>%
    mutate(group = "normal") %>%
    mutate(tree = 1:35)

weights <- rbind(weights_fast, weights_normal)

p <- ggplot(weights, aes(x = tree, y = weight, fill = group)) +
        geom_col(position = "dodge") +
        xlab("Tree") +
        ylab("Weight") +
        theme_bw()

# Save PDF coz ggsave can't handle character spacing properly
savePDF <- function (p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(p)
  invisible(grDevices::dev.off())
}

savePDF(p, here("figures/weights.pdf"), width = 6, height = 4)


