project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.7 Cross Cancer/artifacts/images")
setwd(working_dir)
library(pROC)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source(paste0(project_dir, "./r/utils/ggplot_theme_publication.R"))

set2 <- brewer.pal(3,"Set2")
red <- set2[2]
blue <- set2[3]

metric <- "AUROC"

get_xcomparison_for_model <- function(model_name) {
  xcomparison_results <- list()
  for (train_cancer in cancers) {
    xcomparison <- readRDS(paste0("../", train_cancer, "_cross_comparison_run_data.Rdata"))
    mean_auc <- c()
    sd_auc <- c()
    for (test_cancer in cancers) {
      aucs <- c()
      for (i in 1:10) {
        labels <- xcomparison$labels[[test_cancer]][[i]]
        preds <- xcomparison$preds[[test_cancer]][[model_name]][[i]]
        # aucs <- c(pROC::auc(pROC::roc(labels, preds, direction = "<")), aucs)  
        # aucs <- c(roc.curve(scores.class0 = preds, weights.class0 = labels)$auc, aucs)
        if (metric == "AUROC") {
          aucs <- c(roc.curve(scores.class0 = preds, weights.class0 = labels)$auc, aucs)
        }
        if (metric == "AUPRC") {
          aucs <- c(pr.curve(scores.class0 = preds, weights.class0 = labels)$auc.integral, aucs)
        }
      }
      mean_auc <- c(mean_auc, mean(aucs))
      sd_auc <- c(sd_auc, sd(aucs))
    }
    xcomparison_results[[paste0(train_cancer, "mean")]] <- mean_auc
    xcomparison_results[[paste0(train_cancer, "sd")]] <- sd_auc
  }
  xcomparison_results
}

get_LOCO_results_for_model <- function(model) {
  results <- c()
  for (test_cancer in cancers) {
    data <- readRDS(paste0("../", test_cancer, "_LOCO_run_data.Rdata"))
    aucs <- c()
    for (i in 1:10) {
      labels <- data$labels[[i]]
      preds <- data$preds[[model]][[i]]
      # aucs <- c(pROC::auc(pROC::roc(labels, preds, direction = "<")), aucs)  
      # aucs <- c(roc.curve(scores.class0 = preds, weights.class0 = labels)$auc, aucs)
      aucs <- c(pr.curve(scores.class0 = preds, weights.class0 = labels)$auc.integral, aucs) 
    }
    results[[paste0(test_cancer, "mean")]] <- mean(aucs)
    results[[paste0(test_cancer, "sd")]] <- sd(aucs)
  }
  results
}

# Start plotting
cancers <- c("BRCA", "COAD", "LUAD", "OV")

l0l2 <- get_xcomparison_for_model("L0L2")
l0l2 <- dplyr::bind_rows(l0l2[c(1,3,5,7)])
dimnames(l0l2) <- list(cancers, cancers)

elasticnet <- get_xcomparison_for_model("Elastic Net")
elasticnet <- dplyr::bind_rows(elasticnet[c(1,3,5,7)])
dimnames(elasticnet) <- list(cancers, cancers)

muvr <- get_xcomparison_for_model("MUVR")
muvr <- dplyr::bind_rows(muvr[c(1,3,5,7)])
dimnames(muvr) <- list(cancers, cancers)

rrf <- get_xcomparison_for_model("Random Forest")
rrf <- dplyr::bind_rows(rrf[c(1,3,5,7)])
dimnames(rrf) <- list(cancers, cancers)

xc_d <- dplyr::bind_rows(list(l0l2, elasticnet, muvr, rrf))
xc_d <- stack(xc_d)
xc_d$values <- round(xc_d$values, 2)
xc_d$Tested <- rep(rep(c("BRCA", "COAD", "LUAD", "OV"), 2), 4)
colnames(xc_d) <- c("AUROC", "Trained On", "Tested On")
xc_d$model <- rep(c(rep("L0L2", 4), rep("Elastic Net", 4), rep("MUVR", 4), rep("RRF", 4)), 4)
xc_d$experiment <- "Cross-Cancer"

# get LOCO results
l0l2_loo <- get_LOCO_results_for_model("L0L2")
l0l2_loo <- l0l2_loo[c(1,3,5,7)]
elasticnet_loo <- get_LOCO_results_for_model("Elastic Net")
elasticnet_loo <- elasticnet_loo[c(1,3,5,7)]
muvr_loo <- get_LOCO_results_for_model("MUVR")
muvr_loo <- muvr_loo[c(1,3,5,7)]
rrf_loo <- get_LOCO_results_for_model("Random Forest")
rrf_loo <- rrf_loo[c(1,3,5,7)]

loo_d <- data.frame(round(unname(c(unlist(l0l2_loo), unlist(elasticnet_loo), unlist(muvr_loo), unlist(rrf_loo))), 2), 
                    rep("Others", 8), 
                    rep(c("BRCA", "COAD", "LUAD", "OV"), 4), 
                    c(rep("L0L2", 4), rep("Elastic Net", 4), rep("MUVR", 4), rep("RRF", 4)), 
                    c(rep("LOCO", 16)))
colnames(loo_d) <- c("AUROC", "Trained On", "Tested On", "model", "experiment")


plot_xcomparison_heatmap_2 <- function(x) {
  colnames(x)[2:3] <- c("Trained_on", "Tested_on")
  
  plt <- ggplot(data = x, aes(Tested_on, Trained_on, fill = AUROC))+
    geom_tile(color = "white", height = 1, width=1) +
    geom_text(aes(label=AUROC), size=4) + 
    theme_Publication() +
    theme(
      legend.title = element_blank(),
      legend.key.width=unit(1,"cm"),
      plot.title = element_blank(),
      plot.margin=grid::unit(c(0,0,0,2), "mm"),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()) +
    scale_fill_gradient2(low = blue, high = red, mid = "white", midpoint = 0.5, space = "Lab", limit = c(0,1), 
                         name="AUROC") +
    facet_grid(rows=vars(experiment), cols=vars(model), scales="free_y", space = "free_y")
  return(plt)
}

all_d <- rbind(dplyr::filter(loo_d, grepl('L0L2|MUVR', model)), dplyr::filter(xc_d, grepl('L0L2|MUVR', model)))
all_d$`Trained On` <- factor(all_d$`Trained On`, levels = c("Others", "BRCA", "COAD", "LUAD", "OV"))
all_d$model <- factor(all_d$model, levels = c("L0L2", "MUVR", "L0L2 (LOCO)", "MUVR (LOCO)"))
all_p <- plot_xcomparison_heatmap_2(all_d)
all_p

ggsave(paste0("~/repos/SBSL-modelling-and-analysis/results/figures/heatmaps_mini_", metric, ".pdf"), all_p, width = 16, height = 12, units = "cm")

all_d <- rbind(loo_d, xc_d)
all_d$`Trained On` <- factor(all_d$`Trained On`, levels = c("Others", "BRCA", "COAD", "LUAD", "OV"))
all_d$model <- factor(all_d$model, levels = c("L0L2", "Elastic Net", "MUVR", "RRF", "L0L2 (LOCO)", "Elastic Net (LOCO)", "MUVR (LOCO)", "RRF (LOCO)"))
all_p <- plot_xcomparison_heatmap_2(all_d)
all_p

ggsave(paste0("~/repos/SBSL-modelling-and-analysis/results/figures/heatmaps_all_", metric, ".pdf"), all_p, width = 30, height = 12, units = "cm")




