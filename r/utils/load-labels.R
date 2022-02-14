library(readr)
library(readxl)
library(dplyr)

labels.load <- function(dataset, all_data = FALSE) {
  if (dataset == "isle") {
    labels <- read_excel("~/repos/SBSL-modelling-and-analysis/raw_data/labels/ISLE.xlsx", skip = 2)
    labels <- dplyr::rename(labels, cancer_type = `cancer type tested`)
    labels <- labels[!(labels$`gene1 perturbation` == "mut" | labels$`gene2 perturbation` == "mut"), ]
  } 
  
  if (dataset == "discoversl") {
    labels <- read_delim("~/repos/SBSL-modelling-and-analysis/raw_data/labels/DiscoverSL_trainingSet.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
    labels <- dplyr::rename(labels, cancer_type = Cancer, gene1 = Gene1, gene2 = Gene2, SL = Class)
    labels$SL <- sapply(labels$SL, function(x) ifelse(x == "positive", 1, 0))
  }
  
  if (all_data) {
    return(labels)
  } else {
    return(labels[, c("gene1", "gene2", "cancer_type", "SL")])
  }
  
} 
