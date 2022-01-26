project_dir <- "~/repos/msc-thesis-project/"
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
        aucs <- c(pROC::auc(pROC::roc(labels, preds, direction = "<")), aucs)  
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
      aucs <- c(pROC::auc(pROC::roc(labels, preds, direction = "<")), aucs)  
    }
    results[[paste0(test_cancer, "mean")]] <- mean(aucs)
    results[[paste0(test_cancer, "sd")]] <- sd(aucs)
  }
  results
}





# Start plotting


# get cross cancer results
cancers <- c("BRCA", "COAD", "LUAD", "OV")

l0l2 <- get_xcomparison_for_model("L0L2")
l0l2 <- dplyr::bind_rows(l0l2[c(1,3,5,7)])
dimnames(l0l2) <- list(cancers, cancers)

muvr <- get_xcomparison_for_model("MUVR")
muvr <- dplyr::bind_rows(muvr[c(1,3,5,7)])
dimnames(muvr) <- list(cancers, cancers)

xc_d <- dplyr::bind_rows(list(l0l2, muvr))
xc_d <- stack(xc_d)
xc_d$values <- round(xc_d$values, 2)
xc_d$Tested <- rep(rep(c("BRCA", "COAD", "LUAD", "OV"), 2), 4)
colnames(xc_d) <- c("AUROC", "Trained On", "Tested On")
xc_d$model <- rep(c(rep("L0L2", 4), rep("MUVR", 4)), 4)
xc_d$experiment <- "Cross-Cancer"

# get LOCO results

l0l2_loo <- get_LOCO_results_for_model("L0L2")
l0l2_loo <- l0l2_loo[c(1,3,5,7)]
muvr_loo <- get_LOCO_results_for_model("MUVR")
muvr_loo <- muvr_loo[c(1,3,5,7)]

loo_d <- data.frame(round(unname(c(unlist(l0l2_loo), unlist(muvr_loo))), 2), 
                    rep("Others", 8), 
                    rep(c("BRCA", "COAD", "LUAD", "OV"), 2), 
                    c(rep("L0L2", 4), rep("MUVR", 4)), 
                    c(rep("LOCO", 8)))
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
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()) +
    scale_fill_gradient2(low = blue, high = red, mid = "white", midpoint = 0.5, space = "Lab", limit = c(0,1), 
                         name="AUROC") +
    facet_grid(rows=vars(experiment), cols=vars(model), scales="free_y", space = "free_y")
  return(plt)
}

all_d <- rbind(loo_d, xc_d)
all_d$`Trained On` <- factor(all_d$`Trained On`, levels = c("Others", "BRCA", "COAD", "LUAD", "OV"))
all_d$model <- factor(all_d$model, levels = c("L0L2", "MUVR", "L0L2 (LOCO)", "MUVR (LOCO)"))
all_p <- plot_xcomparison_heatmap_2(all_d)
all_p












print_table <- function(means, sds) {
  means <- round(unlist(means), 2)
  sds <- round(unlist(sds), 2)
  values <- paste0(means, "\u00B1", sds)
  names(values) <- cancers
  values
}


plot_LOCO_heatmap_2 <- function(x_loo, name) {
  melted_cormat_loo <- melt(round(unlist(x_loo), 2))
  ggplot(data = melted_cormat_loo, aes(cancers, 1, fill = value))+
    geom_tile(color = "white", height=1, width=1) +
    geom_text(aes(label=value), size=4) +  
    ggtitle(name) +
    theme_Publication()+ 
    theme(
      legend.title = element_blank(),
      legend.key.width=unit(1,"cm"),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()) +
    coord_fixed() + 
    scale_fill_gradient2(low = blue, high = red, mid = "white", midpoint = 0.5, space = "Lab", limit = c(0,1), 
                         name="AUROC") +
    facet_wrap(~model)
}



logr_loo <- get_LOCO_results_for_model("L0L2")
logr_loo <- logr_loo[c(1,3,5,7)]
logr_loo <- get_LOCO_results_for_model("Elastic Net")
logr_loo <- logr_loo[c(1,3,5,7)]
rf_loo <- get_LOCO_results_for_model("MUVR")
rf_loo <- rf_loo[c(1,3,5,7)]
rf_loo <- get_LOCO_results_for_model("Random Forest")
rf_loo <- rf_loo[c(1,3,5,7)]
p1 <- plot_LOCO_heatmap_2(rf_loo, "MUVR")
p3 <- plot_LOCO_heatmap_2(logr_loo, "Elastic Net")
p2 <- plot_LOCO_heatmap_2(logr_loo, "L0L2")
p4 <- plot_LOCO_heatmap_2(rf_loo, "Random Forest")

sub1 <- ggarrange(p2 + theme(legend.position = "none"), p1 + theme(legend.position = "none"), ncol = 1)
sub1

gg1 <- ggarrange(z2, z1, sub1, common.legend = T, legend = "right", nrow = 1)
ggsave("heatmaps.pdf", gg1, width = 30, height = 11, units = "cm")

xcancergg <- ggarrange(z2, z3, z1, z4, common.legend = T, legend = "bottom", nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
ggsave("Xcancer_all.pdf", xcancergg, width = 30, height = 30, units = "cm")

locogg <- ggarrange(p2, p3, p1, p4, common.legend = T, legend = "bottom", nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
ggsave("LOCO_all.pdf", locogg, width = 30, height = 15, units = "cm")


for (model in c("MUVR", "Random Forest", "Elastic Net", "L0L2")) {
  res <- get_xcomparison_for_model(model)
  res_mean <- dplyr::bind_rows(res[c(1,3,5,7)])
  res_sd <- dplyr::bind_rows(res[c(2,4,6,8)])
  write.table(print_table(res_mean, res_sd), file=paste0(model, "_cross_cancer.tsv"), quote=FALSE, sep='\t', col.names = NA)
}

for (model in c("MUVR", "Random Forest", "Elastic Net", "L0L2")) {
  res <- get_LOCO_results_for_model(model)
  res_mean <- dplyr::bind_rows(res[c(1,3,5,7)])
  res_sd <- dplyr::bind_rows(res[c(2,4,6,8)])
  write.table(print_table(res_mean, res_sd), file=paste0(model, "_LOCO.tsv"), quote=FALSE, sep='\t', col.names = NA)
}




