project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.2 Per Cancer/artifacts/")
setwd(working_dir)
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/monitoring.R"))
source(paste0(project_dir, "r/experiments/baseline/pca-gCMF.R"))
source(paste0(project_dir, "r/experiments/baseline/daisy.R"))
source(paste0(project_dir, "r/experiments/baseline/discoversl.R"))

library(ggplot2)
library(caret)
library(doParallel)
library(foreach)
library(ROCR)
library(pROC)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

print(monitor.system_details())
print(monitor.report_available_workers())
cancer <- args[1]
cancers <- c("BRCA", "COAD", "LUAD", "OV")

for (cancer in cancers) {
  data <- train.get_dataset("combined", cancer)
  preProcValues <- caret::preProcess(data[5:ncol(data)], method = c("center", "scale", "nzv"), uniqueCut=5)
  print(cancer)
  for (v in preProcValues$method$remove) {
    print(paste(v, (length(unique(data[[v]])) - 1), "/", nrow(data) ,(length(unique(data[[v]])) - 1)/nrow(data)))
  }
}


