here::i_am("scripts/plot_weights.R")
library(here)
library(dplyr)
library(magrittr)
library(ggplot2)

# Load weights from .iqtree
file <- "../output/mast_model_6/mast_model_final_6.iqtree"

# Extract weights in one line
weights <- as.numeric(unlist(strsplit(trimws(gsub("Tree weights:|,", "", readLines(file, warn = FALSE)[grep("^Tree weights", readLines(file, warn = FALSE))])), "\\s+")))

weights <- reshape2::melt(weights, value.name = "weight") %>%
    mutate(group = "normal") %>%
    mutate(tree = c(30, 8, 31, 10, 2, 21))

tree = c(30, 8, 31, 10, 2, 21)

p <- ggplot(weights, aes(x = tree, y = weight)) +
  geom_col(position = "dodge") +
  xlab("Tree") +
  ylab("Weight") +
  theme_bw()
p

weights$tree <- factor(weights$tree, levels = weights$tree[order(-weights$weight)])

p <- ggplot(weights, aes(x = tree, y = weight)) +
        geom_col(position = "dodge") +
        xlab("Tree") +
        ylab("Weight") +
        theme_bw()
p

#lowest BIC is the 6-tree model, ie. the model with the 6 trees with the heighest weight
top6_weights <- weights[order(-weights$weight), ]

weight_plot <- ggplot(top6_weights, aes(x = factor(tree, levels = unique(tree)), y = weight))+
  geom_col(fill = "grey45", width = 0.9)+
  geom_text(aes(label = round(weight, digits = 3), y = 0.02),color = "white", size = 5, vjust = 0.5) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2))+
  labs(y = "Model weight", x = "tree")+
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank(),
      #  axis.line.y = element_blank(),
      #  axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 9))+
  guides(colour = 'none')
weight_plot

ggsave(here("figures/weight_6_trees_plot.png"), plot = weight_plot, device = "png", width = 100, height = 50, units = "mm", dpi = 400, bg = "white")


# Save PDF coz ggsave can't handle character spacing properly
savePDF <- function (weight_plot, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(weight_plot)
  invisible(grDevices::dev.off())
}

savePDF(weight_plot, here("figures/weights6tree.pdf"), width = 5, height = 2)


