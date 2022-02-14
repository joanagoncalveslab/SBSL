library(readr)
library(readxl)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/global-vars.R")


# load kegg
kegg_pathway <- read_delim("~/repos/SBSL-modelling-and-analysis/raw_data/msigdb/c2.cp.kegg.v7.0.symbols.gmt", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           col_types = cols(X2 = col_skip()), trim_ws = TRUE)
rownames(kegg_pathway) <- kegg_pathway$X1
kegg_pathway$X1 <- NULL

# load reactome
reactome_pathway <- read_delim("~/repos/SBSL-modelling-and-analysis/raw_data/msigdb/c2.cp.reactome.v7.0.symbols.gmt", 
                               "\t", escape_double = FALSE, col_types = cols(X2 = col_skip()), 
                               trim_ws = TRUE, col_names = FALSE)
rownames(reactome_pathway) <- reactome_pathway$X1
reactome_pathway$X1 <- NULL

# load pid
pid_pathway <- read_delim("~/repos/SBSL-modelling-and-analysis/raw_data/msigdb/c2.cp.pid.v7.0.symbols.gmt", 
                          "\t", escape_double = FALSE, col_types = cols(X2 = col_skip()), 
                          trim_ws = TRUE, col_names = FALSE)
rownames(pid_pathway) <- pid_pathway$X1
pid_pathway$X1 <- NULL

pathway_genes <- na.omit(unique(c(unlist(kegg_pathway), unlist(reactome_pathway), unlist(pid_pathway))))

N <- dim(kegg_pathway)[1] + dim(reactome_pathway)[1] + dim(pid_pathway)[1]
x <- matrix(0, length(pathway_genes), N) == 1
rownames(x) <- pathway_genes

i <- 0
for (g in pathway_genes) {
  in_kegg <- apply(kegg_pathway, 1, function(z) g %in% z)
  in_pid <- apply(pid_pathway, 1, function(z) g %in% z)
  in_reactome <- apply(reactome_pathway, 1, function(z) g %in% z)
  x[g,] <- c(in_kegg, in_pid, in_reactome)
  i <- i + 1
  if (i%%1000 == 0) { print(paste("completed", i, "of", length(pathway_genes))) }
}

x <- x%*%t(x)

calculate_pathway_coparticipation_pvalue <- function(geneA, geneB) {
  out <- tryCatch({
    # total number of successes (geneA in pathway)
    K <- x[geneA, geneA]
    # total number of draws (geneB in pathway)
    n <- x[geneB, geneB]
    # total number of observed sucesses (geneA and geneB has mutation)
    k <- x[geneA, geneB]
    # hypergeometic test
    phyper(k - 1, K, N - K, n, lower.tail = FALSE, log.p = FALSE)
  }, error = function(e) 1
  , warning = function (e) 1)
  
  return(out)
}



for (C in cancer_types) {
  # select genes of interest from labels
  labels <- labels.load(labels_source)
  labels <- dplyr::filter(labels, labels$cancer_type == C)
  unique_genes <- as.matrix(labels[c("gene1", "gene2")])
  unique_genes <- union(unique_genes[,1], unique_genes[,2])
  labelled_pathway_genes <- intersect(unique_genes, pathway_genes)
  
  pathway_values <- c()
  for (i in  1:dim(labels)[1]) {
    geneA <- labels[[i, "gene1"]]
    geneB <- labels[[i, "gene2"]]
    p <- calculate_pathway_coparticipation_pvalue(geneA, geneB)
    pathway_values <- c(pathway_values, p)
    if (i%%200 == 0) { print(paste("completed", i, "of", dim(labels)[1])) }
  }
  
  write.table(pathway_values, paste0(output_dir(), C, "_pathway_coparticipation.txt"), sep = "\t", col.names=c("p.value"), row.names = FALSE)
}
