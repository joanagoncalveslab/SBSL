project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "results")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/ggplot_theme_publication.R"))

cancers <- c("BRCA", "COAD", "LUAD", "OV")
databases = c("ISLE", "discoverSL")

all_genes <- c()
num_labelled_pairs <- 0
for (cancer in cancers) {
  d1 <- train.get_dataset("combined", cancer)[, c(1, 2, 4)]
  all_genes <- c(all_genes, d1$gene1, d1$gene2)
  count_unique_genes <- length(unique(c(d1$gene1, d1$gene2)))
  print(nrow(d1)/choose(count_unique_genes, 2))
  num_labelled_pairs <- nrow(d1) + num_labelled_pairs
}
count_unique_genes <- length(unique(all_genes))
print(count_unique_genes)
print(num_labelled_pairs/(4*choose(count_unique_genes, 2)))


library(ggplot2)
library(gplots)
library(reshape2)
library(ggdendro)
library(RColorBrewer)
library(ComplexHeatmap)
library(dendextend)

convert_to_matrix <- function(d) {
  g <- union(d$gene1, d$gene2)
  m <- matrix(-1, nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(d[1:2])] <- (d$SL == "Y") * 1
  m[as.matrix(d[2:1])] <- (d$SL == "Y") * 1
  
  char_m <- matrix("No Label", nrow = length(g), ncol = length(g));
  colnames(char_m) <- g
  row.names(char_m) <- g
  char_m[as.matrix(d[1:2])] <- ifelse((d$SL == "Y"), "Positive", "Negative")
  char_m[as.matrix(d[2:1])] <- ifelse((d$SL == "Y"), "Positive", "Negative")
  
  return(list(m=m, char_m=char_m))
}

# trying to make nice heatmap
basic_heatmap <- function(cancer, cutoff, show_names=TRUE, height=NA, show_legend=TRUE, special="none") {
  d <- train.get_dataset("combined", cancer)[, c(1, 2, 4)]
  m <- convert_to_matrix(d)
  tbl <- table(c(d$gene1, d$gene2))
  # get row labels
  if (cutoff != 0){
    cutoff_value <- quantile(tbl, cutoff)
    tbl80 <- tbl > cutoff_value
    indexes <- which(rownames(m$m) %in% names(tbl80[tbl80]))
    labels <- rownames(m$m)[indexes]
  } else {
    row_labels <- rownames(m$m)
  }

  paired.cols <- brewer.pal(12,"Paired")
  colors = structure(c("white", paired.cols[6], paired.cols[2]), names = c("-1", "0", "1"))
  
  row_dend <- as.dendrogram(hclust(d = dist(x = m$m)))
  row_order <- order.dendrogram(row_dend)
  row_ordered_labels <- rownames(m$m)[row_order]
  
  ha <- rowAnnotation(Counts = anno_barplot(as.vector(tbl[rownames(m$m)]), 
                                            gp=gpar(border =NA,fill="grey",lty="blank", fontsize=6),
                                            axis_param = list(gp = gpar(fontsize=6))), 
                      annotation_name_gp= gpar(fontsize = 6))

  if (cutoff != 0) {
    ha <- rowAnnotation(Counts = anno_barplot(as.vector(tbl[rownames(m$m)]), 
                                              gp=gpar(border =NA,fill="grey",lty="blank", fontsize=6),
                                              axis_param = list(gp = gpar(fontsize=6))), 
                        Gene = anno_mark(indexes, labels, which="row", labels_gp = gpar(fontsize=6), extend = 1.5), 
                        annotation_name_gp= gpar(fontsize = 6))
  }
  
  # Src Family Kinases: Src, Fyn, Yes, Fgr, Lck, Hck, Blk, Lyn, Frk
  # Non-receptor tyrosine-protein: ABL
  # Receptor tyrosine kinase: EPHA2, RET
  # Signal transducers and activators of transcription (STATs) are activated by Src family kinases in addition to growth factor receptors. STAT5B
  # High affinity nerve growth factor receptor NTRK1
  
  if (special == "mini") {
    highlight_genes <- c("LCK", "ABL2", "FYN", "ABL1", "YES1", "KIT", "EPHA2")
    highlight_row <- rownames(m$m) %in% highlight_genes
    highlight_row_cols <- ifelse(highlight_row, "red", "black")
    fontsize <- 4
  } else {
    highlight_row_cols <- rep("black", nrow(m$m))
    fontsize <- 4
  }
  
  
  
  hmap <- Heatmap(m$m, 
                  col=colors, 
                  height = nrow(m$m)*unit(height, "mm"),
                  show_column_dend = FALSE,
                  cluster_columns = row_dend,
                  cluster_rows = row_dend,
                  row_names_gp = grid::gpar(fontsize = fontsize, col = highlight_row_cols),
                  column_names_gp = grid::gpar(fontsize = 4),
                  show_column_names = FALSE, 
                  show_row_names = show_names,
                  right_annotation = ha,
                  show_heatmap_legend = show_legend,
                  heatmap_legend_param = list(
                    labels = c("No Label", "Negative", "Positive"),
                    title="SL",
                    border="grey",
                    nrow = 3,
                    labels_gp = gpar(fontsize = 6),
                    title_gp = gpar(fontsize = 8),
                    just = c("right", "top")
                  )) 
  
  print(hmap)
  return(list(data=d, matrix=m, hmap=hmap, tbl=tbl, row_dend=row_dend))
}


cancer <- "BRCA"
results <- basic_heatmap(cancer, 0, show_names = TRUE, 2.2, show_legend=FALSE)
pdf(paste0("figures/", cancer, ".dendo-heatmap.pdf"), width = 8, height = 100)
draw(results$hmap)
dev.off()

cancer <- "OV"
results <- basic_heatmap(cancer, 0, show_names = TRUE, 1.25, show_legend=T, special="mini")
pdf(paste0("figures/", cancer, ".dendo-heatmap-mini.pdf"), width = 5.12, height = 5.12)
draw(results$hmap)
dev.off()

cancer <- "OV"
results <- basic_heatmap(cancer, 0, show_names = TRUE, 2.2, show_legend=FALSE)
pdf(paste0("figures/", cancer, ".dendo-heatmap.pdf"), width = 8, height = 8)
draw(results$hmap)
dev.off()

cancer <- "COAD"
results <- basic_heatmap(cancer, 0, show_names = TRUE, 2.2, show_legend=FALSE)
pdf(paste0("figures/", cancer, ".dendo-heatmap.pdf"), width = 8, height = 150)
draw(results$hmap)
dev.off()

cancer <- "LUAD"
results <- basic_heatmap(cancer, 0, show_names = TRUE, 2.2, show_legend=FALSE)
pdf(paste0("figures/", cancer, ".dendo-heatmap.pdf"), width = 8, height = 100)
draw(results$hmap)
dev.off()
# ggsave("OV.dendo-heatmap.pdf", results$hmap, width = 13, height = 9, units = "cm")

paired.cols <- brewer.pal(12,"Paired")
sl_legend  <- Legend(
  labels = c("No Label", "Negative", "Positive"),
  title="SL",
  border="grey",
  ncol = 3,
  labels_gp = gpar(fontsize = 6),
  title_gp = gpar(fontsize = 8),
  legend_gp = gpar(fill = c("white", paired.cols[6], paired.cols[2]))
)
pdf(paste0("figures/sl_legend.pdf"), width = 4, height = 0.5)
draw(sl_legend)
dev.off()

