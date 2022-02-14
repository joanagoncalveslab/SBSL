library(pROC)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(foreach)

library(RColorBrewer)
paired.cols <- brewer.pal(12,"Paired")

source("~/repos/SBSL-modelling-and-analysis/r/utils/ggplot_theme_publication.R")


average_roc_curve_comparisons <- function(d, filename, title) {
  V = 10
  rocs <- list()
  i_aucs <- i_prcs <- list()
  target.fpr <- seq(0, 1, 0.01)
  target.sp <- 1 - target.fpr
  target.sens <- seq(0, 1, 0.01)
  
  for (nm in attributes(d)$names) {

    i_rocs <- foreach(i = 1:V) %do% {roc(d[[nm]][[i]]$labels, d[[nm]][[i]]$predictions, direction = "<")}
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
  
  # plot average ROCs
  ggroc <- dplyr::bind_rows(rocs)
  ggroc$AUC_Model[ggroc$AUPRC_Model == "Random Forest"] <- "RRF"
  ggroc$AUPRC_Model[ggroc$AUPRC_Model == "Random Forest"] <- "RRF"
  ggroc$AUC_Model <- factor(ggroc$AUC_Model, levels = c("L0L2", "Elastic Net", "MUVR", "RRF", "GCATSL", "pca-gCMF", "GRSMF", "DAISY", "DiscoverSL"), ordered = T)
  ggroc$AUPRC_Model <- factor(ggroc$AUPRC_Model, levels = c("L0L2", "Elastic Net", "MUVR", "RRF", "GCATSL", "pca-gCMF", "GRSMF", "DAISY", "DiscoverSL"), ordered = T)
  
  
  rocp <- ggplot(data=ggroc, aes(x=FPR, y=TPR, group=AUC_Model)) +
    geom_line(aes(color=AUC_Model), key_glyph = "rect") + 
    xlim(0, 1) + 
    ylim(0, 1) +
    geom_ribbon(aes(ymin=sapply(TPR-TPR_sdev, function(x) max(0,x)),ymax=sapply(TPR+TPR_sdev, function(x) min(1,x)),fill=AUC_Model), alpha=0.2) +
    ggtitle(paste0(title)) +
    scale_color_manual(values=c(paired.cols[c(1,2,3,4,6,8,10,12,9)], "#000000")) + 
    scale_fill_manual(values=c(paired.cols[c(1,2,3,4,6,8,10,12,9)], "#000000")) + theme_Publication() + theme(
      legend.title = element_blank(),
      aspect.ratio = 1
    )
  
  prcp <- ggplot(data=ggroc, aes(x=Recall, y=Precision, group=AUPRC_Model)) +
    geom_line(aes(color=AUPRC_Model)) + 
    xlim(0, 1) + 
    ylim(0, 1) +
    geom_ribbon(aes(ymin=sapply(Precision-Prec_sdev, function(x) max(0,x)),ymax=sapply(Precision+Prec_sdev, function(x) min(1,x)),fill=AUPRC_Model), alpha=0.2) +
    ggtitle(paste0(title)) +
    scale_color_manual(values=c(paired.cols[c(1,2,3,4,6,8,10,12,9)], "#000000")) + 
    scale_fill_manual(values=c(paired.cols[c(1,2,3,4,6,8,10,12,9)], "#000000")) + theme_Publication() + theme(
      legend.title = element_blank(),
      aspect.ratio = 1
    )
  
  return(list(rocp = rocp, prcp = prcp))
  
}
