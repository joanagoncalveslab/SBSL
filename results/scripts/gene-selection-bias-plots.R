setwd("~/repos/SBSL-modelling-and-analysis/results/")

library(readr)
library(dplyr)
library(xtable)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

paired.cols <- brewer.pal(12,"Paired")
log.cols <- brewer.pal(8,"Accent")
source("../r/utils/ggplot_theme_publication.R")

tested_cancers <- c("BRCA", "LUAD", "COAD", "OV")
singles <- c()
shos <- c()
dhos <- c()
# read all the data
for (cancer in tested_cancers) {
  singles[[cancer]] <- read_csv(paste0(cancer, "_single.csv"))
  shos[[cancer]] <- read_csv(paste0(cancer, "_sho.csv"))
  if (cancer != "OV"){
    dhos[[cancer]] <- read_csv(paste0(cancer, "_dho.csv"))
  }
}


# read all 
ggdatas <- c()
boxps <- c()
for (cancer in tested_cancers){
  if (cancer == "OV") {
    ggdatas[[cancer]] <- data.frame(
      Method = rep(singles[["OV"]]$method[1:70], 2),
      Dropout = factor(c(rep("None", 70), rep("Single", 70)),  levels = c("None", "Single", "Double")),
      AUROC = c(singles[["OV"]]$auroc[1:70], shos[["OV"]]$auroc[1:70]),
      AUPRC = c(singles[["OV"]]$auprc[1:70], shos[["OV"]]$auprc[1:70]),
      AP_3 = c(singles[["OV"]]$ap_3[1:70], shos[["OV"]]$ap_3[1:70]),
      Cancer = "OV"
    )
  } else {
    ggdatas[[cancer]] <- data.frame(
      Method = rep(singles[[cancer]]$method[1:70], 3),
      Dropout = factor(c(rep("None", 70), rep("Single", 70), rep("Double", 70)),  levels = c("None", "Single", "Double")),
      AUROC = c(singles[[cancer]]$auroc[1:70], shos[[cancer]]$auroc[1:70], dhos[[cancer]]$auroc[1:70]),
      AUPRC = c(singles[[cancer]]$auprc[1:70], shos[[cancer]]$auprc[1:70], dhos[[cancer]]$auprc[1:70]),
      AP_3 = c(singles[[cancer]]$ap_3[1:70], shos[[cancer]]$ap_3[1:70], dhos[[cancer]]$ap_3[1:70]),
      Cancer = cancer
    )
  }
  ggdatas[[cancer]]$Method[ggdatas[[cancer]]$Method == "Random Forest"] <- "RRF"
}


data <- rbind(ggdatas[["BRCA"]], ggdatas[["COAD"]], ggdatas[["LUAD"]], ggdatas[["OV"]])
data$Method <- factor(data$Method, levels = c("L0L2", "Elastic Net", "MUVR", "RRF", "GCATSL", "pca-gCMF", "GRSMF"), ordered = T)

plt <- ggplot(data, aes(x=Dropout, y=AUROC, color=Method)) +
  xlab("Dropout") + ylim(0,1) +
  geom_boxplot() + 
  scale_color_manual(values=paired.cols[c(1,2,3,4,6,8,10,12)]) + 
  facet_grid(~Cancer) + theme_Publication() + theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    plot.margin=grid::unit(c(0,0,0,2), "mm"),
  ) + guides(color = guide_legend(nrow = 1))
plt
ggsave("./figures/gene-dropout.pdf", plt, width = 32, height = 9, units = "cm")

plt <- ggplot(data, aes(x=Dropout, y=AUPRC, color=Method)) +
  xlab("Dropout") + ylim(0,1) +
  geom_boxplot() + 
  scale_color_manual(values=paired.cols[c(1,2,3,4,6,8,10,12)]) + 
  facet_grid(~Cancer) + theme_Publication() + theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    plot.margin=grid::unit(c(0,0,0,2), "mm"),
  ) + guides(color = guide_legend(nrow = 1))
plt
ggsave("./figures/gene-dropout-AUPRC.pdf", plt, width = 32, height = 9, units = "cm")

plt <- ggplot(data, aes(x=Dropout, y=AP_3, color=Method)) +
  xlab("Dropout") + ylim(0,1) +
  geom_boxplot() + 
  scale_color_manual(values=paired.cols[c(1,2,3,4,6,8,10,12)]) + 
  facet_grid(~Cancer) + theme_Publication() + theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    plot.margin=grid::unit(c(0,0,0,2), "mm"),
  ) + guides(color = guide_legend(nrow = 1))
plt
ggsave("./figures/gene-dropout-AP_3.pdf", plt, width = 32, height = 9, units = "cm")


project_dir <- "~/repos/SBSL-modelling-and-analysis/"
source(paste0(project_dir, "r/utils/train-model.R"))

melt_data <- function(labels1, d1, labels2, d2, g) {
  m <- matrix("None", nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(labels1[1:2])] <- d1 
  m[as.matrix(labels1[2:1])] <- d1
  m[as.matrix(labels2[1:2])] <- ifelse(m[as.matrix(labels2[1:2])] == d1, "Both", d2)
  m[as.matrix(labels2[2:1])] <- ifelse(m[as.matrix(labels2[2:1])] == d1, "Both", d2)
  m[m == "None"] <- "No Label"
  m <- matrix(factor(m, levels = c("NA", "Both", "DiscoverSL", "Isle")), nrow = length(g))
  melted <- melt(m)
  return(melted)
}

gen_heatmap_2 <- function(melted) {
  red <- paired.cols[6]
  grey <- log.cols[8]
  blue <- paired.cols[2]
  p1 <- ggplot(melted, aes(x = Var2, y = Var1, fill = as.factor(value))) + 
    geom_tile() +
    scale_fill_manual(breaks = c("Both", "DiscoverSL", "Isle"), values =  c(blue, red, grey, "white")) + 
    theme_Publication() + 
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y=element_blank()) + 
    facet_wrap(~Cancer, scales = "free") +
    guides(color = guide_legend(override.aes = list(size = 0.1)))
  return(p1)
}

print("LUAD: Train on DiscoverSL, Test on ISLE")
isle_LUAD <- train.get_dataset("isle", "LUAD")[c(1,2,4)]
discoversl_LUAD <- train.get_dataset("discoversl", "LUAD")[c(1,2,4)]
genes <- Reduce(union, c(isle_LUAD$gene1, isle_LUAD$gene2, discoversl_LUAD$gene1, discoversl_LUAD$gene2))
d1 <- melt_data(discoversl_LUAD, "DiscoverSL", isle_LUAD, "Isle", genes)
d1$Cancer <- "LUAD"


print("BRCA: Train on ISLE, Test on DiscoverSL")
isle_BRCA <- train.get_dataset("isle", "BRCA")[c(1,2,4)]
discoversl_BRCA <- train.get_dataset("discoversl", "BRCA")[c(1,2,4)]
genes <- Reduce(union, c(isle_BRCA$gene1, isle_BRCA$gene2, discoversl_BRCA$gene1, discoversl_BRCA$gene2))
d2 <- melt_data(discoversl_BRCA, "DiscoverSL", isle_BRCA, "Isle", genes)
d2$Cancer <- "BRCA"

melted <- rbind(d1, d2)
hmp <- gen_heatmap_2(melted)
hmp


ggsave("figures/xdata.pdf", hmp, width = 16, height = 11, units = "cm")










