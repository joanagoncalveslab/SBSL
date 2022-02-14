project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.2 Per Cancer/artifacts/images/regularisation")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
library(caret)
library(pROC)
library(foreach)
library(ggpubr)
library(xtable)


elastic_net_regularisation <- function(models, cancer) {
  alphas <- unlist(models$`Elastic Net`[[1]]$results[1])
  lambdas <- unlist(models$`Elastic Net`[[1]]$results[2])
  mean_ROC <- rowMeans(sapply(models$`Elastic Net`, function(x) {unlist(x$results[3])}))
  d <- data.frame(alphas, lambdas, mean_ROC)
  p <- ggplot(d, aes(x=alphas, y=lambdas)) + 
    geom_tile(aes(fill=mean_ROC)) + 
    scale_fill_gradient(low = "white", high = "red") + 
    ggtitle(cancer)
  plot(p)
  ggsave(paste0("EN_", cancer, ".png"), p)
}

RRF_regularisation <- function(models, cancer) {
  mtry <- unlist(models$`Random Forest`[[1]]$results[1])
  coefReg <- unlist(models$`Random Forest`[[1]]$results[2])
  mean_ROC <- rowMeans(sapply(models$`Random Forest`, function(x) {unlist(x$results[3])}))
  d <- data.frame(mtry, coefReg, mean_ROC)
  p <- ggplot(d, aes(x=mtry, y=coefReg)) + 
    geom_tile(aes(fill=mean_ROC)) + 
    scale_fill_gradient(low = "white", high = "red") + 
    ggtitle(cancer)
  plot(p)
  ggsave(paste0("RRF_", cancer, ".png"), p)
}

L0L2_regularisation <- function(models, cancer) {
  gammas <- c()
  lambdas <- c()
  logistic_losses <- c()
  
  for (i in 1:10) {
    g <- models$L0L2[[i]]$fit$gamma
    g <- sapply(g, function(x) rep(x, 50))
    dim(g) <- c(length(g), 1)
    gammas <- c(g, gammas)
    
    cvMeans <- models$L0L2[[i]]$cvMeans
    cvMeans <- sapply(cvMeans, function(x) x)
    dim(cvMeans) <- c(1000, 1)
    logistic_losses <- c(cvMeans, logistic_losses)
    
    
    l <- models$L0L2[[i]]$fit$lambda
    l <- unlist(l)
    lambdas <- c(l, lambdas)
  }
  
  d <- data.frame(gamma = g, lambda = l, logistic_loss = logistic_losses)
  
  p <- ggplot(d, aes(x=gamma, y=lambda, z = logistic_loss)) + 
    stat_summary_2d() + 
    scale_fill_gradient(low = "red", high = "white") + 
    ggtitle(cancer)
  plot(p)
  ggsave(paste0("L0L2_", cancer, ".png"), p)
}

MUVR_feature_selection <- function(models, cancer) {
  nVar <- sapply(models$MUVR, function(x) x$nVar)
  nVarMeans <- rowMeans(nVar)
  nVarSds <- apply(nVar, 1, sd)
  
  aucs <- sapply(models$MUVR, function(x) x$auc[, 2])
  aucMeans <- rowMeans(aucs)
  aucsSds <- apply(aucs, 1, sd)
  
  nVars <- paste0(round(nVarMeans, 2), "+", round(nVarSds, 2))
  aucs <- paste0(round(aucMeans, 2), "+", round(aucsSds, 2))
  
  d <- data.frame(nVars, aucs)
  sink(paste0("MUVR_", cancer, ".txt"))
  xtable(d)
  sink()
}



cancers <- c("BRCA", "COAD", "LUAD", "OV")
plots <- list()
tables <- list()
for (cancer in cancers) {
  predictions <- list()
  for (w in c("single", "baselines")) {
    p <- readRDS(paste0("../../", cancer, "_", w, "_predictions.Rdata"))
    predictions$labels <- p$labels
    predictions$preds <- c(predictions$preds, p$preds)
  }
  print(cancer)
  # predictions$preds[["RRF"]] <- predictions$preds$`Random Forest`
  # predictions$preds$`Random Forest` <- NULL
  for (p in predictions$preds) { print(length(p))}
  
  title <- cancer
  filename <- cancer
  preds <- predictions$preds
  labels <- predictions$labels
  V <- length(preds$`pca-gCMF`)
  
  
  model_names <- names(predictions$preds)[1:4]
  models <- list()
  for (nm in model_names) {
    models[[nm]] <- readRDS(paste0("../../", cancer,"_",nm,".Rdata"))
  }
  elastic_net_regularisation(models, cancer)
  L0L2_regularisation(models, cancer)
  RRF_regularisation(models, cancer) 
  MUVR_feature_selection(models, cancer)
}



