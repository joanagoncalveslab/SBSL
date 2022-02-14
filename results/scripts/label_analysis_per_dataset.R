project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "results")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/ggplot_theme_publication.R"))

library(ggplot2)
library(gplots)
library(reshape2)
library(ggdendro)
library(RColorBrewer)
library(ComplexHeatmap)
library(dendextend)

create_matrix <- function(g, isle_labels, dsl_labels) {
  m <- matrix(-1, nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(isle_labels[1:2])] <- 0 
  m[as.matrix(isle_labels[2:1])] <- 0
  m[as.matrix(dsl_labels[1:2])] <- ifelse(m[as.matrix(dsl_labels[1:2])] == 0, 2, 1)
  m[as.matrix(dsl_labels[2:1])] <- ifelse(m[as.matrix(dsl_labels[2:1])] == 0, 2, 1)
  return(m)
}

basic_heatmap <- function(m, cancer, cutoff, show_names=TRUE, height=NA) {
  d <- train.get_dataset("combined", cancer)[, c(1, 2, 4)]
  tbl <- table(c(d$gene1, d$gene2))
  # get row labels
  if (cutoff != 0){
    cutoff_value <- quantile(tbl, cutoff)
    tbl80 <- tbl > cutoff_value
    indexes <- which(rownames(m) %in% names(tbl80[tbl80]))
    labels <- rownames(m)[indexes]
  } else {
    row_labels <- rownames(m)
  }
  

  paired.cols <- brewer.pal(12,"Paired")
  colors = structure(c("white", paired.cols[6], paired.cols[2], paired.cols[4]), names = c("-1", "0", "1", "2"))
  
  row_dend <- as.dendrogram(hclust(d = dist(x = m)))
  row_order <- order.dendrogram(row_dend)
  row_ordered_labels <- rownames(m)[row_order]
  
  ha <- rowAnnotation(Counts = anno_barplot(as.vector(tbl[rownames(m)]), 
                                            gp=gpar(border =NA,fill="grey",lty="blank", fontsize=3),
                                            axis_param = list(gp = gpar(fontsize=4))), 
                      # Gene = anno_mark(indexes, labels, which="row", labels_gp = gpar(fontsize=4), extend = 1.5),
                      annotation_name_gp= gpar(fontsize = 4))
  
  highlight_row_cols <- rep("black", nrow(m))
  fontsize <- 4
  
  if (any(m==2)) {
    legend_labels <- c("No Label", "ISLE", "DiscoverSL", "Both")
  } else {
    legend_labels <- c("No Label", "ISLE", "DiscoverSL")
  }
  
  hmap <- Heatmap(m, 
                  col=colors, 
                  height = nrow(m)*unit(height, "mm"),
                  show_column_dend = FALSE,
                  cluster_columns = row_dend,
                  cluster_rows = row_dend,
                  row_names_gp = grid::gpar(fontsize = fontsize, col = highlight_row_cols),
                  column_names_gp = grid::gpar(fontsize = 3),
                  show_column_names = FALSE, 
                  show_row_names = TRUE,
                  right_annotation = ha,
                  show_heatmap_legend = FALSE,
                  heatmap_legend_param = list(
                    labels = legend_labels,
                    title="SL",
                    border="grey",
                    nrow = length(legend_labels),
                    labels_gp = gpar(fontsize = 6),
                    title_gp = gpar(fontsize = 8)
                  )) 
  
  print(hmap)
  return(list(data=d, matrix=m, hmap=hmap, tbl=tbl, row_dend=row_dend))
}

print("LUAD: Train on DiscoverSL, Test on ISLE")
isle_C <- train.get_dataset("isle", "LUAD")[c(1,2,4)]
discoversl_C <- train.get_dataset("discoversl", "LUAD")[c(1,2,4)]
genes <- Reduce(union, c(isle_C$gene1, isle_C$gene2, discoversl_C$gene1, discoversl_C$gene2))
common_genes <- intersect(union(isle_C$gene1, isle_C$gene2), union(discoversl_C$gene1, discoversl_C$gene2))
print(paste0(length(common_genes), " genes in common in LUAD"))

mat1 <- create_matrix(genes, isle_C, discoversl_C)
cancer <- "LUAD"
results <- basic_heatmap(mat1, cancer, 0, show_names = FALSE, 2.2)
pdf(paste0("figures/", cancer, ".cross_dataset.dendo-heatmap.pdf"), width = 8, height = 100)
draw(results$hmap)
dev.off()



print("BRCA: Train on ISLE, Test on DiscoverSL")
isle_C <- train.get_dataset("isle", "BRCA")[c(1,2,4)]
discoversl_C <- train.get_dataset("discoversl", "BRCA")[c(1,2,4)]
genes <- Reduce(union, c(isle_C$gene1, isle_C$gene2, discoversl_C$gene1, discoversl_C$gene2))
common_genes <- intersect(union(isle_C$gene1, isle_C$gene2), union(discoversl_C$gene1, discoversl_C$gene2))
print(paste0(length(common_genes), "genes in common in BRCA"))

mat2 <- create_matrix(genes, isle_C, discoversl_C)
cancer <- "BRCA"
results <- basic_heatmap(mat2, cancer, 0, show_names = FALSE, 2.2)
pdf(paste0("figures/", cancer, ".cross_dataset.dendo-heatmap.pdf"), width = 8, height = 100)
draw(results$hmap)
dev.off()

paired.cols <- brewer.pal(12,"Paired")
sl_legend  <- Legend(
  labels = c("No Label", "ISLE", "DiscoverSL", "Both"),
  title="SL",
  border="grey",
  ncol = 4,
  labels_gp = gpar(fontsize = 6),
  title_gp = gpar(fontsize = 8),
  legend_gp = gpar(fill = c("white", paired.cols[6], paired.cols[2], paired.cols[4]))
)
pdf(paste0("figures/dataset_legend.pdf"), width = 4, height = 0.5)
draw(sl_legend)
dev.off()






