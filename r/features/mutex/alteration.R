library(readr)
library(dplyr)
library(discover)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types){
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  unique_genes <- as.matrix(labels[c("gene1", "gene2")])
  unique_genes <- union(unique_genes[,1], unique_genes[,2])
  print("cnv and genes loaded")
  
  ##########################################################
  ## Mutation Matrix Creation
  ##########################################################
  
  cnv <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-all_thresholded.by_genes.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  maf <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  
  data.mut <- as.matrix(cnv[, -(1:3)])
  rownames(data.mut) <- cnv$`Gene Symbol`
  colnames(data.mut) <- substr(colnames(data.mut), 1, 15)
  
  data.amps <- (data.mut >= 2) 
  data.dels <- (data.mut <= -2) 
  
  data.mut <- data.amps | data.dels
  
  mutated_genes <- maf[maf$Variant_Classification != "Silent", c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  mutated_genes$Tumor_Sample_Barcode <- substr(mutated_genes$Tumor_Sample_Barcode, 1, 15)
  mutated_genes <- mutated_genes[mutated_genes$Hugo_Symbol %in% cnv$`Gene Symbol`, ]
  mutated_genes <- mutated_genes[mutated_genes$Tumor_Sample_Barcode %in% colnames(data.mut), ]
  mutated_genes <- mutated_genes[!duplicated(mutated_genes), ]
  
  
  for (i in 1:nrow(mutated_genes)) {
    gene <- mutated_genes$Hugo_Symbol[i]
    sample <- mutated_genes$Tumor_Sample_Barcode[i]
    
    data.mut[gene, sample] <- TRUE
    
    if (i%%1000==0) print(paste(i, "of", nrow(mutated_genes)))
  }
  
  N <- ncol(data.mut)
  x <- data.mut[intersect(unique_genes, rownames(data.mut)), ]
  x <- x%*%t(x)
  
  mutex_values <- numeric(nrow(labels))
  for (i in  1:nrow(labels)) {
    mutex_values[i] <-tryCatch({
      geneA <- labels[[i, "gene1"]]
      geneB <- labels[[i, "gene2"]]
      
      # total number of successes (geneA has mutation)
      K <- x[geneA, geneA]
      # total number of draws (geneB has mutation)
      n <- x[geneB, geneB]
      # total number of observed sucesses (geneA and geneB has mutation)
      k <- x[geneA, geneB]
      
      # hypergeometic test
      1 - phyper(k -1, K, N-K, n, lower.tail = FALSE)
    }, error = function (e) {
      0
    }, warning = function (e) {
      0
    })
  }
  
  
  write.table(mutex_values, paste0(output_dir(), C, "_mutex_alt.txt"),
              sep = "\t", col.names=c("mutex_alt"), row.names = FALSE)
  
  print(paste("finished ", C, " analysis"))
}