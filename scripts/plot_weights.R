here::i_am("scripts/plot_weights.R")
library(here)
library(dplyr)
library(magrittr)
library(ggplot2)

# Load weights from .iqtree
weights <- system("grep '^Tree weights' /maps/projects/racimolab/people/fnl270/GitHub/mast_hiv/output/mast_model_35/mast_model_35.iqtree | sed 's/Tree weights: //; s/,//g'", intern = TRUE)
weights <- as.numeric(strsplit(weights, " +")[[1]])

weights <- reshape2::melt(weights, value.name = "weight") %>%
    mutate(group = "normal") %>%
    mutate(tree = 1:35)

p <- ggplot(weights, aes(x = tree, y = weight)) +
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


