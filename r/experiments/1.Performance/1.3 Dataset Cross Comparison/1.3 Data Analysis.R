project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.3 Dataset Cross Comparison/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/ggplot_theme_publication.R"))

library(ggplot2)
library(reshape2)
library(RColorBrewer)

paired.cols <- brewer.pal(12,"Paired")
log.cols <- brewer.pal(8,"Accent")

gen_heatmap <- function(labels, g, cancer) {
  m <- matrix("No data", nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(labels[1:2])] <- (labels$SL == "Y") * 1
  m[as.matrix(labels[2:1])] <- (labels$SL == "Y") * 1
  
  if (cancer == "OV") {
    m <- m[order(m[1,]),]
  }
  
  melted <- melt(m)
  p1 <- ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + 
    guides(fill = guide_legend(override.aes = 
                                  list(size = .5, 
                                       colour="grey"))) +
    scale_fill_manual(values = c(paired.cols[c(6, 2)], "white"), labels = c("Negative", "Positive", "No Label")) +
    ggtitle(cancer) +
    theme_Publication() + 
    theme(aspect.ratio = 1,
          legend.position = "right",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
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
          axis.ticks.y=element_blank())
  p1
}

gen_heatmap_2 <- function(labels1, d1, labels2, d2, g, colors) {
  m <- matrix("None", nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(labels1[1:2])] <- d1 
  m[as.matrix(labels1[2:1])] <- d1
  m[as.matrix(labels2[1:2])] <- ifelse(m[as.matrix(labels2[1:2])] == d1, "Both", d2)
  m[as.matrix(labels2[2:1])] <- ifelse(m[as.matrix(labels2[2:1])] == d1, "Both", d2)
  m[m == "None"] <- "No Label"
  m <- matrix(factor(m, levels = c("No Label", "Both", "DiscoverSL", "Isle")), nrow = length(g))
  # m <- m[1:400, 1:400]
  
  melted <- melt(m)
  melted$value[1] <- "Both"
  print(head(melted))
  p1 <- ggplot(melted, aes(x = Var2, y = Var1, fill = as.factor(value))) + 
    geom_tile() + guides(fill = guide_legend(override.aes = 
                                               list(size = .5, 
                                                    colour="grey"))) +
    scale_fill_manual(labels = c("Both", "DiscoverSL", "Isle", "No Label"), values =  colors) + 
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
          axis.ticks.y=element_blank())
  p1
}

red <- log.cols[6]
grey <- log.cols[1]
blue <- log.cols[5]

print("LUAD: Train on DiscoverSL, Test on ISLE")
isle_C <- train.get_dataset("isle", "LUAD")[c(1,2,4)]
discoversl_C <- train.get_dataset("discoversl", "LUAD")[c(1,2,4)]
genes <- Reduce(union, c(isle_C$gene1, isle_C$gene2, discoversl_C$gene1, discoversl_C$gene2))
common_genes <- intersect(union(isle_C$gene1, isle_C$gene2), union(discoversl_C$gene1, discoversl_C$gene2))
print(paste0(length(common_genes), " genes in common in LUAD"))
hmp1 <- gen_heatmap_2(discoversl_C, "DiscoverSL", isle_C, "Isle", genes, c(grey, red, blue, "white"))
hmp1
ggsave("xdata_LUAD.pdf", hmp1, width = 11, height = 11, units = "cm")


print("BRCA: Train on ISLE, Test on DiscoverSL")
isle_C <- train.get_dataset("isle", "BRCA")[c(1,2,4)]
discoversl_C <- train.get_dataset("discoversl", "BRCA")[c(1,2,4)]
genes <- Reduce(union, c(isle_C$gene1, isle_C$gene2, discoversl_C$gene1, discoversl_C$gene2))
common_genes <- intersect(union(isle_C$gene1, isle_C$gene2), union(discoversl_C$gene1, discoversl_C$gene2))
print(paste0(length(common_genes), "genes in common in BRCA"))
hmp2 <- gen_heatmap_2(discoversl_C, "DiscoverSL", isle_C, "Isle", genes, c(grey, red,  blue, "white"))
ggsave("xdata_BRCA.pdf", hmp2, width = 11, height = 11, units = "cm")

leg <- get_legend(hmp2)

hmp <- ggarrange(hmp1, hmp2, nrow = 1, labels = c("A", "B"), legend = "bottom", common.legend = TRUE, legend.grob = leg)
ggsave("xdata.pdf", hmp, width = 20, height = 11, units = "cm")




p <- list()
# for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
for (cancer in c("OV")) {
  d <- train.get_dataset("combined", cancer)[c(1,2,4)]
  genes <- union(d$gene1, d$gene2)
  
  # genes <- genes[1:86]
  # d <- d[(d$gene1 %in% genes) & (d$gene2 %in% genes),]
  
  print(length(genes))
  p[[cancer]] <- gen_heatmap(d, genes, cancer) + 
    theme(#legend.position="none", 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          legend.text=element_text(size=11),
          panel.border = element_rect(colour = "grey", fill=NA, size=1)) + 
    ggtitle("")
  ggsave(paste0(cancer, ".label-structure.pdf"), p[[cancer]], width = 12, height = 9, units = "cm")
}
pOV <- p$OV + theme(legend.position="right", legend.direction = "vertical")
pLUAD <- p$LUAD + theme(legend.position="right", legend.direction = "vertical")
library(ggpubr)
gs <- ggarrange(p$LUAD, p$OV, nrow = 1, ncol = 2, labels = c("A", "B"), legend = "bottom", common.legend = T)
ggsave("OV.label-structure.pdf", pOV, width = 13, height = 9, units = "cm")
ggsave("labels.pdf", gs, width = 20, height = 11, units = "cm")

for (cancer in c("BRCA", "LUAD", "COAD", "OV")) {
  print(cancer)
  d <- train.get_dataset("combined", cancer)[c(1,2,4)]
  genes <- union(d$gene1, d$gene2)
  print(length(genes))
  
  all_genes <- c(d$gene1, d$gene2)
  all_genes_c <- table(all_genes) 
  print(length(all_genes_c[all_genes_c > 1])/length(genes))
  
  all_genes_frac <- 100 * all_genes_c/dim(d)[1]
  print(length(all_genes_frac[all_genes_frac > .5]))
  
  total_possible <- choose(length(genes), 2)
  perc_labelled <- 100 * dim(d)[1] /total_possible
  print(perc_labelled)
}



library(ggpubr)
gg1 <- ggarrange(p$BRCA + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$COAD + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$LUAD + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$OV + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                              , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 ncol = 2, 
                 nrow = 2, 
                 legend = "bottom", 
                 common.legend = T, 
                 labels = c("BRCA", "COAD", "LUAD", "OV"),
                 font.label=list(color="black",size=9))
ggsave("gene-selection.png", gg1, width = 20, height = 20, units = "cm")




p <- list()
for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
  d <- train.get_dataset("combined", cancer)[c(1,2,4)]
  genes <- table(c(d$gene1, d$gene2))
  genes <- genes[order(genes, decreasing = TRUE)]
  genes <- as.data.frame(genes)
  p[[cancer]] <- ggplot(genes, aes(Var1, Freq)) + geom_col() + ggtitle(cancer)
}

library(ggpubr)
gg1 <- ggarrange(p$BRCA + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$COAD + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$LUAD + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                                , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 p$OV + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
                              , panel.border = element_rect(colour = "grey", fill=NA, size=1)), 
                 ncol = 2, nrow = 2, legend = "bottom", common.legend = T)

gg1
ggsave("gene-counts.png", gg1, width = 17, height = 17, units = "cm")




