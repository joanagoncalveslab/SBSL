library(RTCGAToolbox)
Sys.setenv(TAR = "/bin/tar")
library(readxl)
library(readr)
library(dplyr)

source("~/repos/msc-thesis-project/r/utils/load-labels.R")
source("~/repos/msc-thesis-project/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types){
  # download GISTIC scores
  # gisticDate <- getFirehoseAnalyzeDates(1)
  # destdir <- "/home/colm/repos/msc-thesis-project/raw_data/firehose/"
  # dataset <- getFirehoseData(dataset=C, gistic2Date=gisticDate, clinic=FALSE, GISTIC = TRUE, destdir = destdir)
  
  # get the thresholded values
  cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", C, "-all_thresholded.by_genes.txt"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  unique_genes <- as.matrix(labels[c("gene1", "gene2")])
  unique_genes <- union(unique_genes[,1], unique_genes[,2])
  print("cnv and genes loaded")
  
  row.names(cnv) <- cnv$`Gene Symbol`
  cnv <- cnv[intersect(unique_genes, cnv$`Gene Symbol`),]
  symbols <- cnv$`Gene Symbol`
  cnv <- cnv[,4:dim(cnv)[2]]
  row.names(cnv) <- symbols
  x1 <- as.matrix(cnv)
  
  all_mutex <- matrix(NA, dim(labels)[1], 4)
  N <- dim(x1)[2]
  
  for (v in c(-2, 2)) {
    
    if (v == -2) {
      x <- (x1 <= v)
    } else {
      x <- (x1 >= v)
    }
    x <- x%*%t(x)
    print("co-copy number altered matrix created")
    
    mutex_values <- matrix(NA, dim(labels)[1], 1)
    for (i in  1:dim(labels)[1]) {
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
    
    if (v == -2) {
      mutex_del <- mutex_values
    } else {
      mutex_amp <- mutex_values
    }
    print("Finished.")
  }
  
  all_mutex[, 1] <- mutex_amp
  all_mutex[, 2] <- mutex_del
  
  maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", C, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  maf <- dplyr::select(maf, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
  maf <- dplyr::filter(maf, Variant_Classification != "Silent")
  maf <- dplyr::filter(maf, Hugo_Symbol %in% unique_genes)
  maf <- dplyr::select(maf, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
  
  # build matrix of genes by samples
  x <- as.matrix(table(maf[c("Hugo_Symbol", "Tumor_Sample_Barcode")])) > 0
  x <- x%*%t(x)
  
  mutex_mut <- matrix(NA, dim(labels)[1], 1)
  for (i in  1:dim(labels)[1]) {
    mutex_mut[i] <- tryCatch({
      geneA <- labels[[i, "gene1"]]
      geneB <- labels[[i, "gene2"]]
      # total number of successes (geneA has mutation)
      K <- x[geneA, geneA]
      # total number of draws (geneB has mutation)
      n <- x[geneB, geneB]
      # total number of observed sucesses (geneA and geneB has mutation)
      k <- x[geneA, geneB]
      
      # hypergeometic test
      1 - phyper(k - 1, K, N-K, n, lower.tail = FALSE)
    }, error = function (e) {
      0
    }, warning = function (e) {
      0
    })
  }
  
  all_mutex[, 3] <- mutex_mut
  
  mutex <- matrix(NA, dim(labels)[1], 1)
  for (i in 1:nrow(mutex)) {
    mutex[i] <- tryCatch({
      metap::sumlog(all_mutex[i, 1:3])$p
    }, error = function(e) {
      1
    }, warning = function(e) {
      1
    })
  }
  
  all_mutex[, 4] <- mutex
  write.table(all_mutex, paste0(output_dir(), C, "_discoverSL_mutex.txt"),
              sep = "\t", col.names=c("mutex_amp", "mutex_del", "mutex_mut", "mutex"), row.names = FALSE)
  
  print(paste("finished ", C, " analysis"))
}