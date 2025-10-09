here::i_am("scripts/plot_ic.R")
library(here)
library(dplyr)
library(magrittr)
library(ggplot2)

OUT <- here("output")
aic <- numeric(35L)
aicc <- numeric(35L)
bic <- numeric(35L)

# This doesnt work for me since i am on windows.. 
# for(i in 2:35){
#     aic[i] <- as.numeric(system(paste0("grep '^Akaike information criterion (AIC) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Akaike information criterion (AIC) score: //; s/,//g'"), intern = TRUE))
#     aicc[i] <- as.numeric(system(paste0("grep '^Corrected Akaike information criterion (AICc) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Corrected Akaike information criterion (AICc) score: //; s/,//g'"), intern = TRUE))
#     bic[i] <- as.numeric(system(paste0("grep '^Bayesian information criterion (BIC) score' ",OUT,"/mast_model_",i,"/mast_model_final_",i,".iqtree | sed 's/Bayesian information criterion (BIC) score: //; s/,//g'"), intern = TRUE))
# }

#using for-loop to get IC values instead. 
for (i in 2:35) {
  # Build the filename
  file <- file.path(OUT, paste0("mast_model_", i), paste0("mast_model_final_", i, ".iqtree"))
  
  # Read all lines
  lines <- readLines(file, warn = FALSE)
  
  # Extract AIC
  aic_line <- grep("^Akaike information criterion \\(AIC\\) score", lines, value = TRUE)
  if (length(aic_line) > 0) {
    aic[i] <- as.numeric(sub("Akaike information criterion \\(AIC\\) score: *|,", "", aic_line))
  }
  
  # Extract AICc
  aicc_line <- grep("^Corrected Akaike information criterion \\(AICc\\) score", lines, value = TRUE)
  if (length(aicc_line) > 0) {
    aicc[i] <- as.numeric(sub("Corrected Akaike information criterion \\(AICc\\) score: *|,", "", aicc_line))
  }
  
  # Extract BIC
  bic_line <- grep("^Bayesian information criterion \\(BIC\\) score", lines, value = TRUE)
  if (length(bic_line) > 0) {
    bic[i] <- as.numeric(sub("Bayesian information criterion \\(BIC\\) score: *|,", "", bic_line))
  }
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

p


p_bic <- ggplot(bic, aes(x = n, y = value)) +
  geom_point(colour = "black", size = 1.2) +
  theme_bw() +
  theme(axis.ticks = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(size = 12, colour = "black"))

p_bic

# Save PDF coz ggsave can't handle character spacing properly
savePDF <- function (p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = "CP1250")
  print(p)
  invisible(grDevices::dev.off())
}

savePDF(p, here("figures/bic.pdf"), width = 4.5, height = 2.8)

write.csv(bic, file = here("output/bic_values.csv"), row.names = FALSE)
