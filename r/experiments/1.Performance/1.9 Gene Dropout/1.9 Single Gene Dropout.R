project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.9 Gene Dropout/artifacts/")
setwd(working_dir)
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/monitoring.R"))
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/experiments/baseline/pca-gCMF.R"))

library(caret)
library(doParallel)
library(foreach)
library(tidyverse)

print(monitor.system_details())
print(monitor.report_available_workers())
nCores <- min(detectCores(), 4)

args <- commandArgs(trailingOnly=TRUE)
cancer <- args[1]
V <- 10

if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 1
  cancer <- "BRCA"
}
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

counts <- models <- preds <- labels <- list()
for (i in 1:V) {
  fold_start <- Sys.time()
  data <- train.balance_cancers(train.get_dataset("combined", cancer), random_seed=Sys.time())
  
  # Build training set where it doesn't have any gene2 from the test set, but all gene1s in the test set appear 
  # in the training set
  common_genes <- intersect(data$gene1, data$gene2)
  test_genes <- data$gene2[!data$gene2 %in% common_genes]
  test <-  data[data$gene2 %in% test_genes, ]
  test_index <- createDataPartition(test$SL, p = (.3*nrow(data))/nrow(test), 
                                     list = FALSE, 
                                     times = 1)
  test <- test[test_index, ]
  # balance training and test sets
  train <- train.balance_cancers(data[!data$gene2 %in% test$gene2, ], random_seed=Sys.time())
  test <- train.balance_cancers(test[test$gene1 %in% train$gene1, ], random_seed=Sys.time())
  
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)
  
  labels[[i]] <- test$SL == "Y"
  
  counts$test_genes[[i]] <- length(unique(test$gene2))
  counts$train_genes[[i]] <- length(unique(train$gene1))
  counts$samples[[i]] <- nrow(train)

  pca_gCMF.preds <- pca_gCMF.predict(train, test, cancer)
  preds[["pca-gCMF"]][[i]] <- pca_gCMF.filter(pca_gCMF.preds, test)
  print("Trained Baselines")
  
  #glmnet
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
  print("Trained Random Forest")
  
  # MUVR
  # https://academic.oup.com/bioinformatics/article/35/6/972/5085367#132378814
  muvr_model <- train.MUVR(train)
  pred <- predict(muvr_model$Fit$rfFitMax, test, type="prob")[,2]
  models[["MUVR"]][[i]] <- muvr_model
  preds[["MUVR"]][[i]] <- pred
  print("Trained MUVR")
  
  print(monitor.get_average_cpu_usage())
  fold_end <- Sys.time()
  print(fold_end - fold_start)
  print(paste0("finished outer fold ", as.character(i)))
}


predictions <- list()
predictions$preds <- preds
predictions$labels <- labels
predictions$counts <- counts
saveRDS(predictions, paste0(cancer, "_predictions_single.Rdata"))

for (nm in names(models)) {
  saveRDS(models[[nm]], paste0(cancer, "_", nm, "_single.Rdata"))
}

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)
