project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.1 Base Datasets/artifacts/")
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

print(monitor.system_details())
print(monitor.report_available_workers())
nCores <- min(detectCores(), 4)
V <- 10
if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 3
}

cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

# variables 
cancers <- train.get_cancer_types()
labels <- preds <- models <- list()

# Repeated stratified nested cross-validation taken from
# https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-10#Sec2

for (i in 1:V) {
  fold_start <- Sys.time()
  data <- train.balance_cancers(train.get_dataset("combined", cancers))
  train_index <- createDataPartition(data$SL, p = .7, 
                                    list = FALSE, 
                                    times = 1)
  train <- data[train_index, ]
  test <- data[-train_index, ]
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)

  # saving results
  labels[[i]] <- test$SL == "Y"
  
  #glmnet
  logr_model <- train.logr(train, f)
  pred <- predict(logr_model, test, type = "prob")[, 2]
  models[["Elastic Net"]][[i]] <- logr_model
  preds[["Elastic Net"]][[i]] <- pred
  
  # L0Learn
  L0L2_model <- train.l0l2(train[5:ncol(train)], train$SL)
  pred <- as.numeric(predict(L0L2_model, newx=data.matrix(test[5:ncol(test)]), lambda=L0L2_model$optimalLambda, gamma=L0L2_model$optimalGamma))
  models[["L0L2"]][[i]] <- L0L2_model
  preds[["L0L2"]][[i]] <- pred
  
  # Random Forest
  rf_model <- train.RRF(train, f)
  pred <- predict(rf_model, test, type = "prob")[, 2]
  models[["Random Forest"]][[i]] <- rf_model
  preds[["Random Forest"]][[i]] <- pred
  
  # MUVR
  # https://academic.oup.com/bioinformatics/article/35/6/972/5085367#132378814
  muvr_model <- train.MUVR(train)
  pred <- predict(muvr_model$Fit$rfFitMax, test, type="prob")[,2]
  models[["MUVR"]][[i]] <- muvr_model
  preds[["MUVR"]][[i]] <- pred

  preds[["DAISY"]][[i]] <- DAISY.predict(test)
  preds[["DiscoverSL"]][[i]] <- discoverSL.predict(test)
  
  print(monitor.get_average_cpu_usage())
  fold_end <- Sys.time()
  print(fold_end - fold_start)
  print(paste0("finished outer fold ", i))

  run_data <- list()
  run_data$preds <- preds
  run_data$models <- models
  run_data$labels <- labels

  saveRDS(run_data, "run_data.Rdata")
}
  
log.dataset_makeup(train, "train")
log.dataset_makeup(test, "test")

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)
