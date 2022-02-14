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


BRCA <- read_csv("BRCA_isle_discoversl_cross_ds.csv")
BRCA$Cancer <- "BRCA"
LUAD <- read_csv("LUAD_discoversl_isle_cross_ds.csv")
LUAD$Cancer  <- "LUAD"

data = rbind(BRCA, LUAD)
data$method[data$method == "Random Forest"] <- "RRF"
names(data)[3] <- "AP_3"
names(data)[1] <- "Method"

data$Method <- factor(data$Method, levels = c("L0L2", "Elastic Net", "MUVR", "RRF", "GCATSL", "pca-gCMF", "GRSMF"), ordered = T)

plt <- ggplot(data, aes(y=ap_3, color=Method)) +
  xlab("Dropout") + ylim(0,1) +
  geom_boxplot() + 
  scale_color_manual(values=paired.cols[c(1,2,3,4,6,8,10,12)]) + scale_x_continuous(breaks = NULL) +
  facet_grid(~Cancer) + theme_Publication() + theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.title = element_blank(),
    plot.margin=grid::unit(c(0,0,0,2), "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) + ggtitle("")
plt
ggsave("./figures/cross-dataset-AP_3", plt, width = 16, height = 10, units = "cm")

