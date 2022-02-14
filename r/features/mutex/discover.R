### https://github.com/NKI-CCB/DISCOVER

library(readr)
library(dplyr)
library(discover)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types) {
  ##########################################################
  ## Mutation Matrix Creation
  ##########################################################
  
  cnv <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-all_thresholded.by_genes.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  maf <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  data.mut <- as.matrix(cnv[, -(1:3)]) != 0
  rownames(data.mut) <- cnv$`Gene Symbol`
  colnames(data.mut) <- substr(colnames(data.mut), 1, 15)
  
  mutated_genes <- maf[maf$Variant_Classification != "Silent", c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  mutated_genes <- mutated_genes[!duplicated(mutated_genes), ]
  mutated_genes$Tumor_Sample_Barcode <- substr(mutated_genes$Tumor_Sample_Barcode, 1, 15)
  mutated_genes <- mutated_genes[mutated_genes$Hugo_Symbol %in% cnv$`Gene Symbol`, ]
  mutated_genes <- mutated_genes[mutated_genes$Tumor_Sample_Barcode %in% colnames(data.mut), ]
  
  for (i in 1:nrow(mutated_genes)) {
    gene <- mutated_genes$Hugo_Symbol[i]
    sample <- mutated_genes$Tumor_Sample_Barcode[i]
    data.mut[gene, sample] <- TRUE
    if (i%%1000==0) print(paste(i, "of", nrow(mutated_genes)))
  }
  
  ##########################################################
  ## Discover Matrix Estimation
  ##########################################################
  events <- discover.matrix(data.mut)
  
  
  labels <- all_labels[all_labels$cancer_type == C, ]
  unique_genes <- as.matrix(labels[c("gene1", "gene2")])
  unique_genes <- union(unique_genes[,1], unique_genes[,2])
  my_subset <- intersect(unique_genes, rownames(events))
  
  result.mutex <- pairwise.discover.test(events[my_subset, ])
  pvalues <- result.mutex$p.values
  
  ##########################################################
  ## Record p-values
  ##########################################################
  mutex <- numeric(nrow(labels)) + 1
  for (j in 1:nrow(labels)) {
    gene1 <- labels$gene1[j]
    gene2 <- labels$gene2[j]
    
    if (gene1 %in% my_subset & gene2 %in% my_subset) {
      out <- tryCatch({
        r <- pvalues[gene1, gene2]
        r <- ifelse(is.na(r), pvalues[gene2, gene1], r)
        ifelse(is.na(r), 1, r)
      }
      , warning = function (e) 1
      , error = function (e) 1)
      
      mutex[j] <- out
    } 
  }
  
  ##########################################################
  ## Write results
  ##########################################################
  write.table(mutex, paste0(output_dir(), C, "_discover_mutex.txt"),
              sep = "\t", col.names=c("discover_mutex"), row.names = FALSE)
  
  print(paste("finished ", C, " analysis"))
}

