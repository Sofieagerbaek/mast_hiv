here::i_am("scripts/plot_bic.R")
library(here)
library(ggplot2)

OUT <- here("output")
bic <- numeric(35L)
for(i in 2:34){
    bic[i] <- as.numeric(system(paste0("grep '^Bayesian information criterion (BIC) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Bayesian information criterion (BIC) score: //; s/,//g'"), intern = TRUE))
}
bic[35] <- as.numeric(system(paste0("grep '^Bayesian information criterion (BIC) score' ",OUT,"/mast_35/mast_35.iqtree | sed 's/Bayesian information criterion (BIC) score: //; s/,//g'"), intern = TRUE))
bic <- bic[-1]
bic <- reshape2::melt(bic, value.name = "BIC")
bic$n <- 2:35

p <- ggplot(bic, aes(x = n, y = BIC)) +
        geom_point() +
        geom_line() +
        theme_bw() +
        xlab("Number of trees") +
        ggtitle("MAST -m GTR+FO+R4+T -fast")

# Save PDF coz ggsave can't handle character spacing properly
savePDF <- function (p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(p)
  invisible(grDevices::dev.off())
}

savePDF(p, here("figures/BIC.pdf"), width = 5, height = 4)
