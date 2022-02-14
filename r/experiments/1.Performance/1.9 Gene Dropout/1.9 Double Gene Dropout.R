project_dir <- "~/repos/SBSL-modelling-and-analysis/"
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
percent_kept <- c(1)
V <- 10

if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 1
  cancer <- "BRCA"
  percent_kept <- c(1)
}
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

counts <- models <- preds <- labels <- list()
for (i in 1:V) {
  fold_start <- Sys.time()
  
  data <- train.balance_cancers(train.get_dataset("combined", cancer), random_seed=Sys.time())
  genes <- union(data$gene1, data$gene2)
  
  # create adjacency matrix for convenience
  adj.genes <- matrix(NA, nrow = length(genes), ncol = length(genes))
  colnames(adj.genes) <- genes
  rownames(adj.genes) <- genes
  for (j in 1:nrow(data)) {
    x <- data[j, ]
    adj.genes[x$gene1, x$gene2] <- adj.genes[x$gene2, x$gene1] <- x$SL
  }
  
  # Get test genes which take up ~ .3 of data set
  # sample training and test genes ensuring a 7:3 train:test ratio when the number of rows for each gene are considered
  train_size <- test_size <- 1
  training_genes <- testing_genes <- c()
  while(length(genes) > 0) {
    if (train_size/test_size <= 4 & length(genes) > 0) {
      training_genes <- c(training_genes, sample(genes, 1))
      genes <- genes[!genes %in% training_genes]
      train_size <- sum(!is.na(adj.genes[training_genes, training_genes]))/2
    }
    if (train_size/test_size > 4 & length(genes) > 0) {
      testing_genes <- c(testing_genes, sample(genes, 1))
      genes <- genes[!genes %in% testing_genes]
      test_size <- sum(!is.na(adj.genes[testing_genes, testing_genes]))/2
    }
  }
  
  train <- data[data$gene1 %in% training_genes & data$gene2 %in% training_genes, ]
  test <-  data[data$gene1 %in% testing_genes & data$gene2 %in% testing_genes, ]

  
  train <- train.balance_cancers(train, random_seed=Sys.time())
  test <- train.balance_cancers(test, random_seed=Sys.time())
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
  f <- train.get_formula(train)
  
  labels[[i]] <- test$SL == "Y"
  
  counts$test_genes[[i]] <- length(unique(test$gene2))
  counts$train_genes[[i]] <- length(unique(train$gene1))
  counts$samples[[i]] <- nrow(train)
  
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

  pca_gCMF.preds <- pca_gCMF.predict(train, test, cancer)
  preds[["pca-gCMF"]][[i]] <- pca_gCMF.filter(pca_gCMF.preds, test)
  print("Trained Baselines")
  
  print(monitor.get_average_cpu_usage())
  fold_end <- Sys.time()
  print(fold_end - fold_start)
  print(paste0("finished outer fold ", as.character(i)))
}


predictions <- list()
predictions$preds <- preds
predictions$labels <- labels
predictions$counts <- counts
saveRDS(predictions, paste0(cancer, "_predictions_double.Rdata"))

for (nm in names(models)) {
  saveRDS(models[[nm]], paste0(cancer, "_", nm, "_double.Rdata"))
}

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)
