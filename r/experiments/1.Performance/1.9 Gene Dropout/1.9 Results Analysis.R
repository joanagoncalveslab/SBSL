project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.9 Gene Dropout/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
library(caret)

cancer <- "BRCA"
data <- readRDS(paste0("../", cancer, "_predictions_single.Rdata"))
title <- cancer
V <- length(data$preds[[1]])
models <- names(data$preds)

i_aucs <- i_prcs <- list()
for (nm in models) {
  i_aucs[[nm]] <- unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = data$preds[[nm]][[i]], weights.class0 = data$labels[[i]])$auc})
  i_prcs[[nm]] <- unlist(foreach(i = 1:V) %do% {pr.curve(scores.class0 = data$preds[[nm]][[i]], weights.class0 = data$labels[[i]])$auc.integral})
}

ggdata_single <- data.frame(
  AUC = unlist(i_aucs), 
  Dropout = "Single",
  Model = unlist(lapply(models, function(x) rep(x, V))))


data <- readRDS(paste0("../", cancer, "_predictions_double.Rdata"))
i_aucs <- i_prcs <- list()
for (nm in models) {
  i_aucs[[nm]] <- unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = data$preds[[nm]][[i]], weights.class0 = data$labels[[i]])$auc})
  i_prcs[[nm]] <- unlist(foreach(i = 1:V) %do% {pr.curve(scores.class0 = data$preds[[nm]][[i]], weights.class0 = data$labels[[i]])$auc.integral})
}
ggdata_double <- data.frame(
  AUC = unlist(i_aucs), 
  Dropout = "Double",
  Model = unlist(lapply(models, function(x) rep(x, V))))

data1 <- readRDS(paste0("../../../1.2 Per Cancer/artifacts/", cancer, "_single_predictions.Rdata"))
data2 <- readRDS(paste0("../../../1.2 Per Cancer/artifacts/", cancer, "_baselines_predictions.Rdata"))
preds <- c(data1$preds, data2$preds)
i_aucs <- i_prcs <- list()
for (nm in models) {
  i_aucs[[nm]] <- unlist(foreach(i = 1:V) %do% {roc.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = data1$labels[[i]])$auc})
  i_prcs[[nm]] <- unlist(foreach(i = 1:V) %do% {pr.curve(scores.class0 = preds[[nm]][[i]], weights.class0 = data1$labels[[i]])$auc.integral})
}
ggdata_none <- data.frame(
  AUC = unlist(i_aucs), 
  Dropout = "None",
  Model = unlist(lapply(models, function(x) rep(x, V))))


ggdata <- dplyr::bind_rows(list(ggdata_single, ggdata_double, ggdata_none))
ggdata$Dropout <- factor(ggdata$Dropout, levels = c("None", "Single", "Double"))

boxp <- ggplot(ggdata, aes(x=Dropout, y=AUC, color=Model)) +
  xlab("Dropout") +
  geom_boxplot() + 
  ggtitle(paste0("AUC score for ", title)) +
  scale_color_manual(values=log.cols) + 
  theme_Publication() + theme(
    legend.title = element_blank()
  )
plot(boxp)

agg_none <- aggregate(AUC ~ Model , ggdata_none, mean)
agg_none_sd <- aggregate(AUC ~ Model , ggdata_none, sd)
agg_single <- aggregate(AUC ~ Model , ggdata_single, mean)
agg_single_sd <- aggregate(AUC ~ Model , ggdata_single, sd)
agg_double <- aggregate(AUC ~ Model , ggdata_double, mean)
agg_double_sd <- aggregate(AUC ~ Model , ggdata_double, sd)

agg <- data.frame(
  Model = agg_none$Model,
  None = paste0(round(agg_none$AUC, 2), "\u00B1", round(agg_none_sd$AUC, 2)),
  Single = paste0(round(agg_single$AUC, 2), "\u00B1", round(agg_single_sd$AUC, 2)),
  Double = paste0(round(agg_double$AUC, 2), "\u00B1", round(agg_double_sd$AUC, 2))
)
write.table(agg, file="gene_dropout.tsv", quote=FALSE, sep='\t', col.names = NA)
