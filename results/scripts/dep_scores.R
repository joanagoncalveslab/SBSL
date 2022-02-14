project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.8 Dep Scores/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
library(doParallel)
library(foreach)
library(caret)
library(tidyverse)
source(paste0(project_dir, "./r/utils/ggplot_theme_publication.R"))


paired.cols <- brewer.pal(8, "Paired")

rocs <- function(preds, labels, cancer) {
  # Threshold averaging of ROC curves as per 
  # https://www.sciencedirect.com/science/article/pii/S016786550500303X?casa_token=2718WOdbeyYAAAAA:4XVSaM_fvwCfa42l6EFOgE4tMsHOCQ9eEq--p4EpMe6aR9IZH0y7SqmfmZuZlB3qpMyqyffbCGE#aep-section-id75
  V <- 10
  rocs <- list()
  i_aucs <- i_prcs <- list()
  target.fpr <- seq(0, 1, 0.01)
  target.sp <- 1 - target.fpr
  target.sens <- seq(0, 1, 0.01)
  for (nm in names(preds)) {
    p <- preds[[nm]]
    
    parts <- str_split(nm, " \\(")[[1]]
    if (length(parts) == 2) {
      model_name <- parts[1]
      experiment <- "No Dep Featrues"
    } else {
      model_name <- parts[1]
      experiment <- "Full Feature Set"
    }
    
    i_rocs <- foreach(i = 1:V) %do% {roc(labels[[i]], preds[[nm]][[i]], direction = "<")}
    i_aucs[[nm]] <- data.frame(
      "AUC" = unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = labels[[i]])$auc}),
      "model" = model_name,
      "experiment" = experiment,
      "cancer" = cancer
    )
  }
  return(i_aucs)
}

format_data <- function(cancer) {
  data <- readRDS(paste0("../", cancer,"_run_data.Rdata"))
  data$preds$RRF <- data$preds$`Random Forest`
  data$preds$`RRF (No Dep)` <- data$preds$`Random Forest (No Dep)`
  data$preds$`Random Forest` <- NULL
  data$preds$`Random Forest (No Dep)` <- NULL
  return(bind_rows(rocs(data$preds, data$labels, cancer)))
}

BRCArocs <- format_data("BRCA")
LUADrocs <- format_data("LUAD")
ggdata <- bind_rows(BRCArocs, LUADrocs)

p <- ggplot(ggdata, aes(x=AUC, y=model, color=model)) +
  ylab("Model") +
  geom_boxplot() + 
  scale_color_manual(values=paired.cols) + xlim(0,1) +
  theme_Publication() +theme(
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.margin=grid::unit(c(0,0,0,2), "mm"),
    legend.position = "bottom"
  ) +
  facet_grid(experiment ~ cancer)
ggsave("~/repos/SBSL-modelling-and-analysis/results/figures/dep_box_plots.pdf", p, width = 16, height = 10, units = "cm")


