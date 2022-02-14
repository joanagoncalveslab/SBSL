project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.5 Embeddings/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/load-raw.R"))

library(tidyr)
library(ggplot2)

cancer <- "BRCA"

dep_scores <- raw.get_ccle_rnai_dependency_scores(cancer)
head(dep_scores)
dim(dep_scores)

# No duplicated genes
genes <- row.names(dep_scores)
duplicated(genes)
any(duplicated(genes))

# counts of which cell lines have the most NAs
na_col_ratio <- colSums(is.na(dep_scores))/nrow(dep_scores)
no_col_df <- data.frame(ratio = na_col_ratio, cell_line = colnames(dep_scores))
ggplot(no_col_df, aes(x = ratio, y = reorder(cell_line, ratio))) +
  geom_col(width = 1) +
  theme(axis.text.y = element_text(size = 4), 
        axis.title.y=element_blank())
ggsave(paste(cancer, "_NA_cell_line_ratio.png"), width = 6, height = 5)
# remove columns with a high number of NAs
no_col_df[no_col_df$ratio > 0.3, ]

# plot little histograms for the distributions of each cell line
ggplot(gather(dep_scores[1:12]), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')
ggsave(paste(cancer, "_cell_line_hist.png"), width = 6, height = 5)

# counts of which genes have NA values
no_genes_na <- rowSums(is.na(dep_scores))/ncol(dep_scores)
no_genes_df <- data.frame(
  ratio = no_genes_na,
  gene = row.names(dep_scores))
ggplot(no_genes_df, aes(x=ratio)) + geom_histogram()
ggsave(paste(cancer, "_NA_gene_ratio.png"), width = 6, height = 5)
# remove genes with mostly NAs
no_genes_df[no_genes_df$ratio > 0.3, ]

ggplot(gather(as.data.frame(t(dep_scores)[,1:12])), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')
ggsave(paste(cancer, "_gene_hist.png"), width = 6, height = 5)

# Analyse for near zero variance
preProcValues <- caret::preProcess(dep_scores, method = c("nzv"))

# plot correlations
cormat <- round(cor(data.matrix(na.omit(dep_scores))), 2)

library(ggplot2)
library(reshape2)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# Melt the correlation matrix
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) +
  coord_fixed()
ggsave("correlation_heatmap.png", width = 8, height = 8, device = "png")

highly_cor <- which(abs(cormat) > .7, arr.ind = T)
highly_cor[apply(highly_cor, 1, function(x) x[1] != x[2]), ]

# analyse for multicollinearity problem
analyse_multicollinearity <- function() {
  summary(glm.fit)
  vifs <- car::vif(glm.fit)
  vifs[vifs > 5]
}

# How to deal with NAs?
# Remove cell lines with very high number of NAs
# For genes with only a few NAs, maybe look at imputing values with Nearest Neighbour? Values are highly correlated
# So this approach might work well
no_na_dep_scores <- na.omit(dep_scores)
dim(no_na_dep_scores)

