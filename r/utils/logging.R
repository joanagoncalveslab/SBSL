library(ROCR)
library(pROC)
library(PRROC)
library(ggplot2)
library(foreach)
library(dplyr)
library(tidyverse)
library(caret)
library(ggpubr)
library(ggrepel)
library(iml)
library(reshape2)
library(RColorBrewer)
library(Metrics)

source(paste0(project_dir, "r/utils/ggplot_theme_publication.R"))


rie <- function (x, y, decreasing = TRUE, alpha = 20) 
{
  if (length(x) != length(y)) {
    stop(paste("The length of scores should be equal to number of labels."))
  }
  N <- length(y)
  n <- length(which(y == 1))
  ord <- order(x, decreasing = decreasing)
  m_rank <- which(y[ord] == 1)
  s <- sum(exp(-alpha * m_rank/N))
  ra <- n/N
  ri <- (N - n)/N
  random_sum <- (n/N) * (1 - exp(-alpha))/(exp(alpha/N) - 1)
  return(s/random_sum)
}


average_precision <- function(x, y) {
  ap <- 0
  n <- length(y)
  K <- floor(n/4)
  ord <- order(x, decreasing = T)
  y <- y[ord]
  for (k in 1:K) {
    ap <- ap + Metrics::precision(y[1:k], numeric(k) + 1)
  }
  ap/K
}

log.cols <- brewer.pal(8,"Set2")

paired.cols <- brewer.pal(8, "Paired")

log.many_cols <- c("#ac4a51",
          "#60c352",
          "#795bd0",
          "#589831",
          "#c943a9",
          "#acb742",
          "#cd74de",
          "#3c9754",
          "#d9407f",
          "#6dc38b",
          "#d13a40",
          "#3fc1bf",
          "#d86632",
          "#5b83e2",
          "#d59d35",
          "#8950a2",
          "#617e33",
          "#b296db",
          "#837a33",
          "#5667a8",
          "#d3a26c",
          "#53a2d5",
          "#9a5c2a",
          "#37825e",
          "#ad5b94",
          "#e18079",
          "#964267",
          "#e386b3")

log.experiment_run <- function () {
  fileConn<-file("last-run.info")
  writeLines(c(
    format(Sys.time(), "%a %b %d %X %Y"),
    Sys.info()["nodename"]
  ), fileConn)
  close(fileConn)
}

log.dataset_makeup <- function(dataset, filename = NULL) {
  dataset_details <- table(cancer = dataset$cancer_type, SL = dataset$SL)
  print(dataset_details)
  write.table(dataset_details, file = paste0(filename, "_details.txt"), sep = "\t", row.names = T)
}

log.logr.regularisation_curves <- function(models, filename, title) {
  results <- list()
  tmp <- models[[1]]$results[c(1,2,3)]
  lambda <- unique(tmp$lambda)
  alpha <- unique(tmp$alpha)
  for (j in 1:length(alpha)) results[[j]] <- list()
  for (i in 1:V) {
    ROCs <- models[[i]]$results[c(1,2,3)]
    for (j in 1:length(alpha)) {
      results[[j]][[i]] <- ROCs$ROC[ROCs$alpha == alpha[j]]
    }
  }
  for (j in 1:length(alpha)) {
    results[[j]][["mean"]] <- apply(do.call(rbind, results[[j]]), 2, mean)
    results[[j]][["sd"]] <- apply(do.call(rbind, results[[j]]), 2, sd)
  }
  
  ggparams <- data.frame(
    ROC_mean = unlist(lapply(results, function(x) x$mean)),
    ROC_sd = unlist(lapply(results, function(x) x$sd)),
    lambda = rep(lambda, length(alpha)),
    alpha = as.factor(unlist(lapply(alpha, function(x) rep(round(x, 2), length(lambda)))))
  )
  
  p <- ggplot(data=ggparams, aes(x=lambda, y=ROC_mean, group=alpha)) +
    geom_line(aes(color=alpha)) + 
    ggtitle(paste0("Effect of Elastic Net regularization parameters for ", title)) + 
    theme(
      plot.title = element_text(size=10)
    ) +
    guides(color = guide_legend(ncol=2)) + 
    labs(color = "Alpha")
  plot(p)
  ggsave(paste0(filename, "_logr_regularisation.png"), width = 6, height = 4)
}

log.rf.regularisation_curves <- function(models, filename, title) {
  mtry <- list()
  mtry_roc <- list()
  for (i in 1:V) {
    tmp <- models[[i]]$results
    mtry <- tmp$mtry
    mtry_roc[[i]] <- tmp$ROC
  }
  ROC_mean <- apply(do.call(rbind, mtry_roc), 2, mean)
  ROC_sd <- apply(do.call(rbind, mtry_roc), 2, sd)
  
  ggparams <- data.frame(
    ROC_mean,
    ROC_sd,
    mtry
  )
  
  p <- ggplot(data=ggparams, aes(x=mtry, y=ROC_mean)) +
    geom_line() + 
    geom_ribbon(aes(ymin=ROC_mean-ROC_sd,ymax=ROC_mean+ROC_sd), alpha=0.2) + 
    ggtitle(paste0("Effect of Random Forest mtry regularization parameter for ", title)) + 
    theme(
      plot.title = element_text(size=10),
      legend.title = element_blank()
    )
  plot(p)
  ggsave(paste0(filename, "_rf_regularisation.png"), width = 6, height = 4)
}

log.print_results <- function(preds, labels, filename, title) {
  # Threshold averaging of ROC curves as per 
  # https://www.sciencedirect.com/science/article/pii/S016786550500303X?casa_token=2718WOdbeyYAAAAA:4XVSaM_fvwCfa42l6EFOgE4tMsHOCQ9eEq--p4EpMe6aR9IZH0y7SqmfmZuZlB3qpMyqyffbCGE#aep-section-id75
  rocs <- list()
  i_aucs <- i_prcs <- i_ncgd <- list()
  target.fpr <- seq(0, 1, 0.01)
  target.sp <- 1 - target.fpr
  target.sens <- seq(0, 1, 0.01)
  for (nm in names(preds)) {
    p <- preds[[nm]]
    i_rocs <- foreach(i = 1:V) %do% {roc(labels[[i]], preds[[nm]][[i]], direction = "<")}
    i_aucs[[nm]] <- unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = labels[[i]])$auc})
    i_prcs[[nm]] <- unlist(foreach(i = 1:V) %do% {pr.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = labels[[i]])$auc.integral})
    i_ncgd[[nm]] <- unlist(foreach(i = 1:V) %do% {average_precision(preds[[nm]][[i]], labels[[i]] * 1)})
  }
  
  # plot AUCs
  ggdata <- data.frame(
    Model = names(preds),
    PRC_mean = sapply(i_prcs, mean),
    PRC_sd = sapply(i_prcs, sd),
    ROC_mean = sapply(i_aucs, mean),
    ROC_sd = sapply(i_aucs, sd),
    nCGD_mean = sapply(i_ncgd, mean),
    nCGD_sd = sapply(i_ncgd, sd)
  )
  ggdata
}

log.average_roc_curve_comparisons <- function(preds, labels, filename, title) {
  # Threshold averaging of ROC curves as per 
  # https://www.sciencedirect.com/science/article/pii/S016786550500303X?casa_token=2718WOdbeyYAAAAA:4XVSaM_fvwCfa42l6EFOgE4tMsHOCQ9eEq--p4EpMe6aR9IZH0y7SqmfmZuZlB3qpMyqyffbCGE#aep-section-id75
  rocs <- list()
  i_aucs <- i_prcs <- list()
  target.fpr <- seq(0, 1, 0.01)
  target.sp <- 1 - target.fpr
  target.sens <- seq(0, 1, 0.01)
  for (nm in names(preds)) {
    p <- preds[[nm]]
    i_rocs <- foreach(i = 1:V) %do% {roc(labels[[i]], preds[[nm]][[i]], direction = "<")}
    i_aucs[[nm]] <- unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = labels[[i]])$auc})
    i_prcs[[nm]] <- unlist(foreach(i = 1:V) %do% {pr.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = labels[[i]])$auc.integral})
    i_tpr <- foreach(i = 1:V) %do% {
      pROC::coords(i_rocs[[i]], target.sp, input = "specificity", ret = c("tpr")) %>% map_df(rev)
    }
    i_precision <- foreach(i = 1:V) %do% {
      pROC::coords(i_rocs[[i]], target.sens, input = "sensitivity", ret = c("precision"))
    }
    tpr <- fpr <- recall <- precision <- list()
    for (i in 1:V) {
      tpr[[i]] <- c(i_tpr[[i]]$tpr, 0)
      fpr[[i]] <- c(target.sp, 0)
      recall[[i]] <- c(target.sens, 0)
      precision[[i]] <- c(i_precision[[i]]$precision, 1)
    } 
    TPR <- apply(dplyr::bind_cols(tpr), 1, mean)
    TPR_sdev <- apply(dplyr::bind_cols(tpr), 1, sd)
    FPR <- apply(dplyr::bind_cols(fpr), 1, mean)
    Precision <- apply(dplyr::bind_cols(precision), 1, mean)
    Precision[is.nan(Precision)] <- 1
    Prec_sdev <- apply(dplyr::bind_cols(precision), 1, sd)
    Prec_sdev[is.nan(Prec_sdev)] <- 0
    Recall <- apply(dplyr::bind_cols(recall), 1, mean)
    rocs[[nm]] <- data.frame(
      TPR, 
      TPR_sdev,
      FPR, 
      Precision,
      Prec_sdev, 
      Recall,
      Model = nm,
      # AUC_Model = paste0(nm, " (", round(mean(i_aucs[[nm]]),2), "±", round(sd(i_aucs[[nm]]),2), ")"),
      # AUPRC_Model = paste0(nm, " (", round(mean(i_prcs[[nm]]),2), "±", round(sd(i_prcs[[nm]]),2), ")")
      AUC_Model = paste0(nm),
      AUPRC_Model = paste0(nm)
      )
  }
  
  # plot AUCs
  ggdata <- data.frame(
    AUC = unlist(i_aucs), 
    Model = unlist(lapply(names(preds), function(x) rep(x, V))))
  boxp <- ggplot(ggdata, aes(x=AUC, y=Model, color=Model)) +
    ylab("Model") +
    geom_boxplot() + 
    ggtitle(paste0("Average AUC score for ", title)) +
    scale_color_manual(values=log.cols) + 
    theme(
      plot.title = element_text(size=9),
      legend.title = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_blank()
    )
  
  # plot average ROCs
  ggroc <- dplyr::bind_rows(rocs)
  rocp <- ggplot(data=ggroc, aes(x=FPR, y=TPR, group=AUC_Model)) +
    geom_line(aes(color=AUC_Model), key_glyph = "rect") + 
    xlim(0, 1) + 
    ylim(0, 1) +
    geom_ribbon(aes(ymin=sapply(TPR-TPR_sdev, function(x) max(0,x)),ymax=sapply(TPR+TPR_sdev, function(x) min(1,x)),fill=AUC_Model), alpha=0.2) +
    ggtitle(paste0(title)) +
    scale_color_manual(values=log.cols) + 
    scale_fill_manual(values=log.cols) + theme_Publication() + theme(
      legend.title = element_blank(),
      aspect.ratio = 1
    )
  
  prcp <- ggplot(data=ggroc, aes(x=Recall, y=Precision, group=AUPRC_Model)) +
    geom_line(aes(color=AUPRC_Model)) + 
    xlim(0, 1) + 
    ylim(0, 1) +
    geom_ribbon(aes(ymin=sapply(Precision-Prec_sdev, function(x) max(0,x)),ymax=sapply(Precision+Prec_sdev, function(x) min(1,x)),fill=AUPRC_Model), alpha=0.2) +
    ggtitle(paste0(title)) +
    scale_color_manual(values=log.cols) + 
    scale_fill_manual(values=log.cols) + theme_Publication() + theme(
      legend.title = element_blank(),
      aspect.ratio = 1
    )
  
  return(list(boxp = boxp, rocp = rocp, prcp = prcp))
}

# Feature importance using IML R Package
pred_fun <- function(model, newdata) {
  predict(model, newdata, type = "prob")
}

pred_fun_L0L2 <- function(model, newdata) {
  as.numeric(predict(model$fit, newx=data.matrix(newdata),
                     lambda=model$optimalLambda,
                     gamma=model$optimalGamma))
}

loss_fun <- function(actual, predicted) {
  1 - Metrics::auc(actual, predicted)
}

calc_imp_over_n_models <- function(train_fun, train, test, f) {
  X <- test[5:ncol(test)]
  y = as.numeric(test$SL == "Y" )
  
  features <- labels(terms(f))[order(labels(terms(f)))]
  imp <- NULL
  for (i in 1:10) {
    m.logr <- train_fun(train, f)
    predictor.logr <- Predictor$new(
      model = m.logr, 
      data = X, 
      y = y, 
      predict.function = pred_fun,
      class = 2,
      type = "prob"
    )
    results <- FeatureImp$new(predictor.logr, loss = loss_fun, n.repetitions = 50)$results
    results <- results[order(results$feature), 2]
    imp <- cbind(imp, results)
  }
  sdev <- apply(imp, 1, sd)
  importance <- apply(imp, 1, mean)
  imp <- data.frame(features, importance, sdev)
  imp <- imp[order(-imp$importance),]
  imp
}

log.plot_var_importance_and_ale_plots <- function(m, test, cancer, L0L2) {
  X <- test[5:ncol(test)]
  y = as.numeric(test$SL == "Y")
  
  predictor <- Predictor$new(
    model = m, 
    data = X, 
    y = y, 
    predict.function = ifelse(L0L2, pred_fun_L0L2, pred_fun),
    type = "prob"
  )
  imp <- FeatureImp$new(predictor, loss = loss_fun, n.repetitions = 50)$results
  imp
  


  
  
  # # Comparison of new variables
  # imp.logr2 <- imp.logr[order(imp.logr$feature), ]
  # imp.rf2 <- imp.rf[order(imp.rf$feature), ]
  # data2 <- data.frame(variables = as.factor(imp.logr2$feature), logr = imp.logr2$importance, rf = imp.rf2$importance)
  # row.names(data2) <- data2$variables
  # p3 <- ggplot(data2[grepl("mutex", data2$variables), ], aes(x=logr, y=rf)) + 
  #   geom_point(color = log.many_cols[1], size = 3) + 
  #   geom_text_repel(aes(label = variables), color = log.many_cols[1], size = 3) +
  #   theme(legend.position = "none") + 
  #   geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed")
  # 
  # p4 <- ggplot(data2[grepl("(RNAi|CRISPR)", data2$variables), ], aes(x=logr, y=rf)) + 
  #   geom_point(color = log.many_cols[2], size = 3) + 
  #   geom_text_repel(aes(label = variables), color = log.many_cols[2], size = 3) +
  #   theme(legend.position = "none") + 
  #   geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  # 
  # p5 <- ggplot(data2[grepl("(gtex|tumour_corr|normal_corr|diff_exp)", data2$variables), ], aes(x=logr, y=rf)) + 
  #   geom_point(color = log.many_cols[3], size = 3) + 
  #   geom_text_repel(aes(label = variables), color = log.many_cols[3], size = 3) +
  #   theme(legend.position = "none") + 
  #   geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  # 
  # p6 <- ggplot(data2[grepl("(logrank)", data2$variables), ], aes(x=logr, y=rf)) + 
  #   geom_point(color = log.many_cols[7], size = 3) + 
  #   geom_text_repel(aes(label = variables), color = log.many_cols[7], size = 3) +
  #   theme(legend.position = "none") + 
  #   geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  # 
  # gg2 <- ggarrange(p3, p4, p5, p6,
  #                  labels = c("C", "D", "E", "F"),
  #                  ncol = 4, nrow = 1)
  # plot(gg2)
  # ggsave(paste0(cancer, "_new_var_comparison.png"), width = 12, height = 4, device = "png")
  # 
  # # Effect Size vs P-Value
  # pvalue_data <- data2[grepl("(gtex|tumour_corr|normal_corr|RNAi|CRISPR|diff_exp)", data2$variables), ]
  # pvalue_data$group <- sapply(pvalue_data$variables, function (x) ifelse(grepl("pvalue", x), "P-Value", "Effect Size"))
  # scatterPlot <- ggplot(pvalue_data,aes(x=logr, y=rf, color=group)) + 
  #   geom_point() + geom_rug() +
  #   scale_color_manual(values = log.many_cols[10:11]) + 
  #   theme(legend.position=c(0.05,0.99), legend.justification=c(0,1))
  # 
  # plot(scatterPlot)
  # ggsave(paste0(cancer, "_effect_vs_pval.png"), width = 4, height = 4, device = "png")
  # 
  # ggarrange(gg1, gg2, ncol = 1, nrow = 2)
  # ggsave(paste0(cancer, "_all_var_imp.png"), width = 9, height = 6, device = "png")
  # 
  # ale.rf <- FeatureEffects$new(predictor.rf)
  # ale.rf$plot()
  # ggsave(paste0(cancer, "_rf_feature_effects.png"), width = 12, height = 9, device = "png")
  # ale.logr <- FeatureEffects$new(predictor.logr)
  # ale.logr$plot()
  # ggsave(paste0(cancer, "_logr_feature_effects.png"), width = 12, height = 9, device = "png")
}

log.dep_score_interactions <- function(models, test, cancer) {
  X <- test[5:ncol(test)]
  y = as.numeric(test$SL == "Y")
  
  predictor.rf <- Predictor$new(
    model = models$MUVR$Fit$rfFitMax, 
    data = X, 
    y = y, 
    predict.function = pred_fun,
    class = 2,
    type = "prob"
  )
  
  plot(Interaction$new(predictor.rf, feature = "CRISPR_dep_stat"))
  ggsave(paste0(cancer, "_rf_interactions_CRISPR_dep_stat.png"), width = 6, height = 4, device = "png")
  plot(Interaction$new(predictor.rf, feature = "RNAi_dep_stat"))
  ggsave(paste0(cancer, "_rf_interactions_RNAi_dep_stat.png"), width = 6, height = 4, device = "png")
}

log.pca_gCMF_stability <- function(preds, cancer) {
  pca_gCMF_aucs <- do.call("rbind", foreach (i = 1:V) %do% {
    foreach (j = 1:10, .combine = "c") %do% {
      pROC::auc(roc(labels[[i]], preds$`pca-gCMF`[[i]][[j]], direction = "<"))
    }
  })
  colnames(pca_gCMF_aucs) <- paste("Run", 1:ncol(pca_gCMF_aucs), sep = "_")
  row.names(pca_gCMF_aucs) <- paste("Fold", 1:nrow(pca_gCMF_aucs), sep = "_")
  pca_gCMF_aucs_melted <- melt(pca_gCMF_aucs)
  colnames(pca_gCMF_aucs_melted) <- c("Fold", "Run", "AUC")
  ggplot(pca_gCMF_aucs_melted, aes(x = Run, y = AUC, group = Fold)) + 
    scale_color_manual(values=log.many_cols) + 
    geom_line(aes(color = Fold)) + 
    geom_hline(yintercept=0.5, color="grey55", linetype="dashed") + 
  ggsave(paste0(cancer, "_pca_gCMF_stability.png"), width = 6, height = 4, device = "png")
  apply(pca_gCMF_aucs, 1, sd)
}
