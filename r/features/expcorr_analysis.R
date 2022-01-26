library(RTCGAToolbox)
Sys.setenv(TAR = "/bin/tar")
library(readxl)
library(readr)
library(dplyr)

source("~/repos/msc-thesis-project/r/utils/load-labels.R")
source("~/repos/msc-thesis-project/r/utils/load-raw.R")
source("~/repos/msc-thesis-project/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types){
  print(paste("starting ", C, " analysis"))
  # download expression data
  # runDate <- getFirehoseRunningDates(1)
  # destdir <- "/home/colm/repos/msc-thesis-project/raw_data/firehose"
  # brcaData <- getFirehoseData(dataset=C, runDate=runDate, clinic=FALSE, RNASeq2GeneNorm = TRUE)
  
  # select genes of interest from labels
  
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  unique_genes <- as.matrix(labels[c("gene1", "gene2")])
  unique_genes <- union(unique_genes[,1], unique_genes[,2])
  
  rna <- raw.get_tcga_expression(C, unique_genes)
  normal_samples <- substr(colnames(rna), 14, 15) == "11"
  
  tumour_corr <- numeric(nrow(labels))
  tumour_corr.pvalue <- numeric(nrow(labels))
  normal_corr <- numeric(nrow(labels))
  normal_corr.pvalue <- numeric(nrow(labels))
  
  for (i in  1:nrow(labels)) {
    t.corr <- tryCatch({
      gene1_exp_tumour <- as.numeric(rna[labels[[i, "gene1"]], !normal_samples])
      gene2_exp_tumour <- as.numeric(rna[labels[[i, "gene2"]], !normal_samples])
      cor.test(gene1_exp_tumour, gene2_exp_tumour, method = "pearson")
    }, error = function(e) {
      c(estimate = 0, p.value = 1)
    }, warning = function(e) {
      c(estimate = 0, p.value = 1)
    })
    
    n.corr <- tryCatch({
      gene1_exp_normal <- as.numeric(rna[labels[[i, "gene1"]], normal_samples])
      gene2_exp_normal <- as.numeric(rna[labels[[i, "gene2"]], normal_samples])
      cor.test(gene1_exp_normal, gene2_exp_normal, method = "pearson")
    }, error = function(e) {
      c(estimate = 0, p.value = 1)
    }, warning = function(e) {
      c(estimate = 0, p.value = 1)
    })
  
    tumour_corr[i] <- t.corr[["estimate"]]
    tumour_corr.pvalue[i] <- t.corr[["p.value"]]
    normal_corr[i] <- n.corr[["estimate"]]
    normal_corr.pvalue[i] <- n.corr[["p.value"]]
    
    if (i%%100 == 0) print(paste("completed", i, "of", nrow(labels)))
  }
  
  corr_values <- data.frame(tumour_corr, tumour_corr.pvalue, normal_corr, normal_corr.pvalue)
  write.table(corr_values, paste0(output_dir(), C, "_exp_corr.txt"), sep = "\t", col.names=colnames(corr_values), row.names = FALSE)
  print(paste("finished ", C, " analysis"))
}
