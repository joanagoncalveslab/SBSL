project_dir <- "~/repos/SBSL-modelling-and-analysis/"
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
nCores <- min(detectCores(), 4)
V <- 10
what_to_train_options <- c("baselines", "unbalanced", "all", "single")
what_to_train <- args[2]

if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 1
  cancer <- "COAD"
  what_to_train <- "all"
}
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

print(paste("Runnning for", cancer, what_to_train))
labels <- preds <- models <- list()
predictions <- list()
# Repeated stratified nested cross-validation taken from
# https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-6-10#Sec2

try(what_to_train %in% what_to_train_options)

train_test_indices <- readRDS("train_test_indices.Rdata")
for (i in 1:V) {
  fold_start <- Sys.time()
  data <- train.get_dataset("combined", cancer)
  train_index <- train_test_indices[[cancer]][[i]]
  train <- data[train_index, ]
  test <- data[-train_index, ]
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  saveRDS(preProcValues, paste0(cancer, "_preProcValues.Rdata"))
  break
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)
    
  test <- train.balance_cancers(test)

  # saving results
  labels[[i]] <- test$SL == "Y"

  if (what_to_train == "baselines") {
    train <- train.balance_cancers(train)
    print(table(train$SL, train$cancer_type))
    preds[["DAISY"]][[i]] <- DAISY.predict(test)
    preds[["DiscoverSL"]][[i]] <- discoverSL.predict(test)
    
    pca_gCMF.preds <- pca_gCMF.predict(train, test, cancer)
    preds[["pca-gCMF"]][[i]] <- pca_gCMF.filter(pca_gCMF.preds, test)
    print("Trained Baselines")
  } 
  
  if (what_to_train == "single") {
    #glmnet
    train <- train.balance_cancers(train)
    print(table(train$SL, train$cancer_type))
    logr_model <- train.logr(train, f)
    pred <- predict(logr_model, test, type = "prob")[, 2]
    models[["Elastic Net"]][[i]] <- logr_model
    preds[["Elastic Net"]][[i]] <- pred
    print("Trained Elastic Net")
    
    # L0Learn
    L0L2_model <- train.l0l2(train[5:ncol(train)], train$SL)
    pred <- as.numeric(predict(L0L2_model, newx=data.matrix(test[5:ncol(test)]), lambda=L0L2_model$optimalLambda, gamma=L0L2_model$optimalGamma))
    models[["L0L2"]][[i]] <- L0L2_model
    preds[["L0L2"]][[i]] <- pred
    print("Trained L0L2")
    
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
    print("Trained MUVR")
  } 
  
  #######################################
  # Unbalanced (All) Cancers
  #######################################
  train <- data[train_index, ]
  test <- data[-train_index, ]
  other_cancers <- cancers[cancers != cancer]
  other_train <- train.get_dataset("combined", other_cancers)
  train <- dplyr::bind_rows(other_train, train)
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)
  test <- train.balance_cancers(test)
  
  if (what_to_train == "unbalanced") {
    train <- train.balance_classes(train)
    print(table(train$SL, train$cancer_type))
    #glmnet
    logr_model <- train.logr(train, f)
    pred <- predict(logr_model, test, type = "prob")[, 2]
    models[["Elastic Net (Unbalanced)"]][[i]] <- logr_model
    preds[["Elastic Net (Unbalanced)"]][[i]] <- pred
    print("Trained Elastic Net")
    
    # L0Learn
    L0L2_model <- train.l0l2(train[5:ncol(train)], train$SL)
    pred <- as.numeric(predict(L0L2_model, newx=data.matrix(test[5:ncol(test)]), lambda=L0L2_model$optimalLambda, gamma=L0L2_model$optimalGamma))
    models[["L0L2 (Unbalanced)"]][[i]] <- L0L2_model
    preds[["L0L2 (Unbalanced)"]][[i]] <- pred
    print("Trained L0L2 (Unbalanced)")
    
    # Random Forest
    rf_model <- train.RRF(train, f)
    pred <- predict(rf_model, test, type = "prob")[, 2]
    models[["Random Forest (Unbalanced)"]][[i]] <- rf_model
    preds[["Random Forest (Unbalanced)"]][[i]] <- pred
    
    muvr_model <- train.MUVR(train)
    pred <- predict(muvr_model$Fit$rfFitMax, test, type="prob")[,2]
    models[["MUVR (Unbalanced)"]][[i]] <- muvr_model
    preds[["MUVR (Unbalanced)"]][[i]] <- pred
    print("Trained MUVR (Unbalanced)")
  } 
  
  if (what_to_train == "all") {
    #######################################
    # Balanced Cancers
    #######################################
    train <- train.balance_cancers(train)
    print(table(train$SL, train$cancer_type))
    #glmnet
    logr_model <- train.logr(train, f)
    pred <- predict(logr_model, test, type = "prob")[, 2]
    models[["Elastic Net (All)"]][[i]] <- logr_model
    preds[["Elastic Net (All)"]][[i]] <- pred
    print("Trained Elastic Net")
    
    # L0Learn
    L0L2_model <- train.l0l2(train[5:ncol(train)], train$SL)
    pred <- as.numeric(predict(L0L2_model, newx=data.matrix(test[5:ncol(test)]), lambda=L0L2_model$optimalLambda, gamma=L0L2_model$optimalGamma))
    models[["L0L2 (All)"]][[i]] <- L0L2_model
    preds[["L0L2 (All)"]][[i]] <- pred
    print("Trained L0L2 (All)")
    
    # Random Forest
    rf_model <- train.RRF(train, f)
    pred <- predict(rf_model, test, type = "prob")[, 2]
    models[["Random Forest (All)"]][[i]] <- rf_model
    preds[["Random Forest (All)"]][[i]] <- pred
    
    muvr_model <- train.MUVR(train)
    pred <- predict(muvr_model$Fit$rfFitMax, test, type="prob")[,2]
    models[["MUVR (All)"]][[i]] <- muvr_model
    preds[["MUVR (All)"]][[i]] <- pred
    print("Trained MUVR (All)")
  }
  
  print(monitor.get_average_cpu_usage())
  fold_end <- Sys.time()
  print(fold_end - fold_start)
  print(paste("finished outer fold ", i, cancer, what_to_train))

  predictions$preds <- preds
  predictions$labels <- labels
  saveRDS(predictions, paste0(cancer, "_", what_to_train ,"_predictions.Rdata"))

  for (nm in names(models)) {
    saveRDS(models[[nm]], paste0(cancer, "_", nm, ".Rdata"))
  }
}

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)


# train_counts <- c()
# test_counts <- c()
# for (i in 1:V) {
#   fold_start <- Sys.time()
#   train_index <- train_test_indices[[cancer]][[i]]
#   train <- data[train_index, ]
#   test <- data[-train_index, ]
#   
#   train <- train.balance_cancers(train)
#   test <- train.balance_cancers(test)
#   
#   train_counts <- c(length(union(train$gene1, train$gene2)), train_counts)
#   test_counts <- c(length(union(test$gene1, test$gene2)), test_counts)
#   
#   print(nrow(train))
#   print(nrow(test))
# }
