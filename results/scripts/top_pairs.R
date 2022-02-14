setwd("~/repos/SBSL-modelling-and-analysis/results/")
library(readr)
library(dplyr)
library(xtable)
library(rjson)
library(ggpubr)
source("./scripts/plot_curves.R")


top5predictions = c()
  
for (cancer in c("BRCA", "LUAD")) {
  predictions <- list()
  p <- fromJSON(paste(readLines(paste0("single_", cancer, "_colm_paper.json")), collapse=""))
  
  for (fold in 1:10) {
    predictions <- p$L0L2[[fold]]$predictions
    po <- order(predictions, decreasing = TRUE)
    p$L0L2[[fold]]$predictions[po][1:5]
    top5predictions = c(p$L0L2[[fold]]$test_names[po][1:5], top5predictions)
    # print(p$L0L2[[fold]]$labels[po][1:5])
  }
  print(cancer)
  print(sort(table(top5predictions), decreasing = TRUE)[1:5])
}


# BRCA (Top 3)
# TP53|BOP1 
# - https://www.hindawi.com/journals/jo/2021/3603030/ (shows BOP1 involved in BRCA and p53 pathway in BRCA)
# - https://journals.asm.org/doi/10.1128/MCB.21.13.4246-4255.2001

# TP53|CRYGS
# No useful papers found

# TP53|EEF1D
# No useful papers found

# LUAD (Top 3)
# KRAS|ANAPC4

# KRAS|FBL

# KRAS|PSMA1