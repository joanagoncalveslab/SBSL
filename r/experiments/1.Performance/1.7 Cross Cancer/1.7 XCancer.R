project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.7 Cross Cancer/artifacts/")
setwd(working_dir)
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/monitoring.R"))
source(paste0(project_dir, "r/experiments/baseline/pca-gCMF.R"))
source(paste0(project_dir, "r/experiments/baseline/daisy.R"))
source(paste0(project_dir, "r/experiments/baseline/discoversl.R"))

# load parrallel processing
library(ggplot2)
library(caret)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)

print(monitor.system_details())
print(monitor.report_available_workers())
train_cancer <- args[1]
nCores <- min(detectCores(), 4)
V <- 10

if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 1
  train_cancer <- "BRCA"
}
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

cancers <- c("BRCA", "COAD", "LUAD", "OV")
dataset <- "combined"

####
# Cross Cancer Comparisons
####
labels <- preds <- models <- list()
for (i in 1:V) {
  fold_start <- Sys.time()
  data <- train.balance_cancers(train.get_dataset("combined", train_cancer), random_seed = NULL)
  train_index <- createDataPartition(data$SL, p = .8, 
                                     list = FALSE, 
                                     times = 1)
  train <- data[train_index, ]
  test <- data[-train_index, ]
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)
  
  models[["Elastic Net"]][[i]] <- train.logr(train, f)
  models[["L0L2"]][[i]] <- train.l0l2(train[5:ncol(train)], train$SL)
  models[["MUVR"]][[i]] <- train.MUVR(train)
  models[["Random Forest"]][[i]] <- train.RRF(train, f)
  
  labels[[train_cancer]][[i]] <- test$SL == "Y"
  preds[[train_cancer]][["Elastic Net"]][[i]] <- predict(models[["Elastic Net"]][[i]], test, type = "prob")[, 2]
  preds[[train_cancer]][["L0L2"]][[i]] <- as.numeric(predict(models[["L0L2"]][[i]], newx=data.matrix(test[5:ncol(test)]),
                                                            lambda=models[["L0L2"]][[i]]$optimalLambda, 
                                                            gamma=models[["L0L2"]][[i]]$optimalGamma))
  preds[[train_cancer]][["MUVR"]][[i]] <- predict(models[["MUVR"]][[i]]$Fit$rfFitMax, test, type="prob")[,2]
  preds[[train_cancer]][["Random Forest"]][[i]] <- predict(models[["Random Forest"]][[i]], test, type = "prob")[, 2]

  for (test_cancer in cancers[cancers != train_cancer]) {
    test <- train.balance_cancers(train.get_dataset(dataset, test_cancer), random_seed = NULL)
    test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
    labels[[test_cancer]][[i]] <- test$SL == "Y"
    preds[[test_cancer]][["Elastic Net"]][[i]] <- predict(models[["Elastic Net"]][[i]], test, type = "prob")[, 2]
    preds[[test_cancer]][["L0L2"]][[i]] <- as.numeric(predict(models[["L0L2"]][[i]], newx=data.matrix(test[5:ncol(test)]),
                                                              lambda=models[["L0L2"]][[i]]$optimalLambda, 
                                                              gamma=models[["L0L2"]][[i]]$optimalGamma))
    preds[[test_cancer]][["MUVR"]][[i]] <- predict(models[["MUVR"]][[i]]$Fit$rfFitMax, test, type="prob")[,2]
    preds[[test_cancer]][["Random Forest"]][[i]] <- predict(models[["Random Forest"]][[i]], test, type = "prob")[, 2]
  }
  
  print(monitor.get_average_cpu_usage())
  fold_end <- Sys.time()
  print(fold_end - fold_start)
  print(paste0("finished outer fold ", i))
}

cross_comparison_run_data <- list()
cross_comparison_run_data$preds <- preds
cross_comparison_run_data$labels <- labels
saveRDS(cross_comparison_run_data, paste0(train_cancer, "_cross_comparison_run_data.Rdata"))

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)