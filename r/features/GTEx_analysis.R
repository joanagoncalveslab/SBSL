library(readr)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/load-raw.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types) {
  ##################################################
  ## Get Ensembl Ids for genes of interest
  ##################################################
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  genes <- as.matrix(labels[c("gene1", "gene2")])
  genes <- union(genes[,1], genes[,2])

  gtex_expression <- raw.get_GTEx_expression(C, genes)
  
  ##################################################
  ## Get gene expression correlation
  ##################################################
  gtex_corr <- numeric(nrow(labels))
  gtex_corr.pvalue <- numeric(nrow(labels))
  for (i in 1:nrow(labels)) {
    gene1 <- labels$gene1[i]
    gene2 <- labels$gene2[i]
    
    r <- tryCatch({
      gene1_exp <- as.numeric(gtex_expression[gene1, ])
      gene2_exp <- as.numeric(gtex_expression[gene2, ])
      cor.test(gene1_exp, gene2_exp, method = "pearson")
    }, error = function(e) {
      c(estimate = 0, p.value = 1)
    }, warning = function(e) {
      c(estimate = 0, p.value = 1)
    })
    
    gtex_corr[i] <- r[["estimate"]]
    gtex_corr.pvalue[i] <- r[["p.value"]]
    
    if (i%%100 == 0) print(paste("completed", i, "of", nrow(labels)))
  }
  gtex_results <- data.frame(gtex_corr, gtex_corr.pvalue)
  write.table(gtex_results, paste0(output_dir(), C, "_gtex.txt"), sep = "\t", col.names=colnames(gtex_results), row.names = FALSE)
  print(paste("finished ", C, " analysis"))
}
