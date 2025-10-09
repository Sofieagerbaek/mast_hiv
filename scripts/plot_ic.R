here::i_am("scripts/plot_ic.R")
library(here)
library(dplyr)
library(magrittr)
library(ggplot2)

OUT <- here("output")
aic <- numeric(35L)
aicc <- numeric(35L)
bic <- numeric(35L)
for(i in 2:35){
    aic[i] <- as.numeric(system(paste0("grep '^Akaike information criterion (AIC) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Akaike information criterion (AIC) score: //; s/,//g'"), intern = TRUE))
    aicc[i] <- as.numeric(system(paste0("grep '^Corrected Akaike information criterion (AICc) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Corrected Akaike information criterion (AICc) score: //; s/,//g'"), intern = TRUE))
    bic[i] <- as.numeric(system(paste0("grep '^Bayesian information criterion (BIC) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Bayesian information criterion (BIC) score: //; s/,//g'"), intern = TRUE))
}
aic <- aic[-1]
aicc <- aicc[-1]
bic <- bic[-1]
aic <- reshape2::melt(aic, value.name = "value") %>% mutate(n = 2:35, estimator = "AIC")
aicc <- reshape2::melt(aicc, value.name = "value") %>% mutate(n = 2:35, estimator = "AICc")
bic <- reshape2::melt(bic, value.name = "value") %>% mutate(n = 2:35, estimator = "BIC")
ic <- rbind(aic, aicc, bic)

p <- ggplot(ic, aes(x = n, y = value, colour = estimator)) +
        geom_point() +
        geom_line() +
        theme_bw() +
        xlab("Number of trees") +
        ggtitle("Model fit")

# Save PDF coz ggsave can't handle character spacing properly
savePDF <- function (p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(p)
  invisible(grDevices::dev.off())
}

savePDF(p, here("figures/ic.pdf"), width = 5, height = 4)
