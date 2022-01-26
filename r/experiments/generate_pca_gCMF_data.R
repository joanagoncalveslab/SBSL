project_dir <- "~/repos/msc-thesis-project/"
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/experiments/baseline/pca-gCMF.R"))

library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)
cancer <- args[1]
nCores <- min(detectCores(), 4)
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)
print(paste0("Runnning for ", cancer))
pca_gCMF.generate_feature_matrices(cancer)