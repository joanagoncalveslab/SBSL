project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.5 Embeddings/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/analyse-model.R"))
source(paste0(project_dir, "r/utils/monitoring.R"))
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/load-raw.R"))
source("../autoencoder.R")

library(ggplot2)
library(caret)
library(doParallel)
library(foreach)
library(pROC)

load_dependency_scores <- function(cancers) {
  dep_scores <- raw.get_ccle_rnai_dependency_scores(cancers)
  dep_scores$gene <- row.names(dep_scores)
  dep_scores <- na.omit(dep_scores)
  row.names(dep_scores) <- dep_scores$gene
  dep_scores$gene <- NULL
  dep_scores
}

encode_data <- function(d, encoder, scores) {
  d <- d[, 1:4]
  gene1.dep_scores <- scores[d$gene1, ]
  gene1.dep_scores_na_i <- unique(which(is.na(gene1.dep_scores), arr.ind = T)[,1])
  gene2.dep_scores <- scores[d$gene2, ]
  gene2.dep_scores_na_i <- unique(which(is.na(gene2.dep_scores), arr.ind = T)[,1])
  na_i <- c(gene1.dep_scores_na_i, gene2.dep_scores_na_i)
  
  d <- d[-na_i, ]
  gene1.dep_scores <- gene1.dep_scores[-na_i, ]
  gene2.dep_scores <- gene2.dep_scores[-na_i, ]
  
  gene1.encoded <- encoder %>% predict(data.matrix(gene1.dep_scores))
  gene2.encoded <- encoder %>% predict(data.matrix(gene2.dep_scores))
  d <- dplyr::bind_cols(list(d, as.data.frame(gene1.encoded), as.data.frame(gene2.encoded)))
  d
}

use_raw_features <- function(d, dep_scores) {
  d <- d[, 1:4]
  gene1.dep_scores <- dep_scores[d$gene1, ]
  gene1.dep_scores_na_i <- unique(which(is.na(gene1.dep_scores), arr.ind = T)[,1])
  gene2.dep_scores <- dep_scores[d$gene2, ]
  gene2.dep_scores_na_i <- unique(which(is.na(gene2.dep_scores), arr.ind = T)[,1])
  na_i <- c(gene1.dep_scores_na_i, gene2.dep_scores_na_i)
  d <- d[-na_i, ]
  gene1.dep_scores <- gene1.dep_scores[-na_i, ]
  gene2.dep_scores <- gene2.dep_scores[-na_i, ]
  d1 <- dplyr::bind_cols(list(data.frame(SL = d$SL), as.data.frame(gene1.dep_scores), as.data.frame(gene2.dep_scores)))
  d2 <- dplyr::bind_cols(list(data.frame(SL = d$SL), as.data.frame(gene2.dep_scores), as.data.frame(gene1.dep_scores)))
  dplyr::bind_rows(list(d1, d2))
}

train_on_embedding_with_n_features <- function(dep_scores, train, test, n) {
  # train autoencoder
  n <- floor(nrow(dep_scores) * .8)
  X.train <- data.matrix(dep_scores[1:n, -ncol(dep_scores)])
  X.test <- data.matrix(dep_scores[n:nrow(dep_scores), -ncol(dep_scores)])
  encoder <- autoencode(X.train, X.test, 20)
  test <- encode_data(test, encoder, dep_scores[-ncol(dep_scores)])
  train <- encode_data(train, encoder, dep_scores[-ncol(dep_scores)])
  logr.m <- train.logr(train[4:ncol(train)], as.factor(SL) ~ .)
}

args <- commandArgs(trailingOnly=TRUE)

print(monitor.system_details())
print(monitor.report_available_workers())
cancer <- args[1]
nCores <- min(detectCores(), 4)
V <- 10

if (Sys.info()["nodename"] == "colm-Inspiron-5577") {
  V <- 10
  cancer <- "BRCA"
}

cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
  scores <- load_dependency_scores(cancer)
  labels <- preds <- models <- list()
  for (i in 1:V) {
    data <- train.balance_cancers(train.get_dataset("combined", cancer))
    train_index <- createDataPartition(data$SL, p = .8, 
                                       list = FALSE, 
                                       times = 1)
    train <- data[train_index, ]
    test <- data[-train_index, ]
    train <- use_raw_features(train, scores)
    test <- use_raw_features(test, scores)
    # data should be centered and scaled for use in L0L2 algorithm 
    preProcValues <- caret::preProcess(train[2:ncol(train)], method = c("center", "scale"))
    train <- cbind(SL = train$SL, predict(preProcValues, train[2:ncol(train)]))
    test <- cbind(test$SL, predict(preProcValues, test[2:ncol(test)]))
    colnames(train)[1] <- colnames(test)[1] <- "SL"
    # generalised linear model fit (for VIF measurements)
    glm.fit <- glm(SL ~ ., data = train, family = "binomial")
    glm.results <- predict(glm.fit, test)
    # L0 L2 penalised logisic regression models
    l0l2.cvfit <- L0Learn.cvfit(data.matrix(train[2:ncol(train)]), train$SL,
                                nFolds=10,
                                loss="Logistic",
                                nGamma = 5,
                                penalty="L0L2")
    cvMeans <- unlist(lapply(l0l2.cvfit$cvMeans, min))
    optimalGammaIndex = which.min(cvMeans)
    optimalLambdaIndex = which.min(l0l2.cvfit$cvMeans[[optimalGammaIndex]])
    l0l2.cvfit$optimalLambda = l0l2.cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
    l0l2.cvfit$optimalGamma = l0l2.cvfit$fit$gamma[optimalGammaIndex]
    l0l2.results <- as.numeric(predict(l0l2.cvfit,
                                       newx=data.matrix(test[2:ncol(test)]),
                                       lambda=l0l2.cvfit$optimalLambda,
                                       gamma=l0l2.cvfit$optimalGamma))
    
    labels[[i]] <- test$SL
    
  
    models[["glm"]][[i]] <- glm.fit
    models[["L0L2"]][[i]] <- l0l2.cvfit
    
    preds[["glm"]][[i]] <- glm.results
    preds[["L0L2"]][[i]] <- l0l2.results
    
    print(monitor.get_average_cpu_usage())
  }
  run_data <- list()
  run_data$preds <- preds
  run_data$models <- models
  run_data$labels <- labels
  saveRDS(run_data, paste0(cancer, "_raw_features_run_data.Rdata"))
}

stopCluster(cl)
log.experiment_run()
registerDoSEQ()
end_time <- Sys.time()
print(end_time - start_time)

# attempt 1. co-dependency matrix embeddings
# codep <- cor(t(dep_scores[-ncol(dep_scores)]))
# codep <- (abs(codep) > 0.5) * 1
# row.names(codep) <- colnames(codep) <- dep_scores$gene
# n <- floor(nrow(codep) * .8)
# X.train <- data.matrix(codep[1:n, ])
# X.test <- data.matrix(codep[n:nrow(codep), ])
# source("../autoencoder.R")
# encoder <- binary_autoencoder(X.train, X.test, 100)
# train <- encode_data(train, encoder, as.data.frame(codep))
# test <- encode_data(test, encoder, data.frame(codep))
# 
# logr.m <- train.logr(train[4:ncol(train)], as.factor(SL) ~ .)
# results <- train.analyse(logr.m, test[4:ncol(test)])
# 
# plot(results$roc)


# attempt 2. embed each gene in a lower dimensional space, and stack vectors







