library(dplyr)
library(readr)
library(survival)
library(survminer)

cancer_types <- c("BRCA", "LUAD", "OV", "COAD", "KIRC", "LAML")

for (C in cancer_types) {
  # runDate <- getFirehoseRunningDates(1)
  # destdir <- "/home/colm/repos/SBSL-modelling-and-analysis/raw_data/firehose"
  # brcaData <- getFirehoseData(dataset=C, runDate=runDate, clinic=TRUE, destdir = destdir)
  
  ##########################################################
  ## Clinical 
  ##########################################################
  clinical <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-",C, "-Clinical.txt"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  cols <- unname(unlist(clinical[, 1]))
  clinical[, 1] <- NULL
  clinical <- as.data.frame(t(clinical))
  colnames(clinical) <- cols
  
  
  ##########################################################
  ## Copy Number Variant 
  ##########################################################
  cnv <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-all_thresholded.by_genes.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  cnv_genes <- unlist(cnv[, 1])
  cnv <- cnv[, -c(1, 2, 3)]
  cnv_samples <- lapply(colnames(cnv), function (x) tolower(substr(x, 1, 12)))
  cnv <- as.data.frame(t(cnv))
  rownames(cnv) <- cnv_samples
  colnames(cnv) <- cnv_genes
  
  
  ##########################################################
  ## Mutation Data 
  ##########################################################
  maf <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  maf <- dplyr::filter(maf, Variant_Classification != "Silent")
  maf <- dplyr::select(maf, Hugo_Symbol, Tumor_Sample_Barcode)
  maf$Tumor_Sample_Barcode <- unlist(lapply(maf$Tumor_Sample_Barcode, function (x) tolower(substr(x, 1, 12))))
  
  
  ##########################################################
  ## Find alterations
  ##########################################################
  alt <- (cnv <= -2) | (cnv >= 2)
  for (i in 1:nrow(maf)) {
    gene <- maf$Hugo_Symbol[i]
    sample <- maf$Tumor_Sample_Barcode[i]
    tryCatch({
      alt[sample, gene] <- alt[sample, gene] | TRUE
    }, error = function(x){})
    if (i%%100 == 0) print(paste(i, "of", nrow(maf)))
  }
  
  ##########################################################
  ## mRNA Data 
  ##########################################################
  rna <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", C, "-RNAseq2GeneNorm.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  rna <- rna[-(1:30), ]
  rna_samples <- unname(sapply(colnames(rna)[-1], function (x) tolower(substr(x, 1, 12))))
  rna_genes <- unname(sapply(unlist(rna[, 1]), function (x) strsplit(x, "\\|")[[1]][1]))
  rna <- rna[, -1]
  rna <- rna[, !duplicated(rna_samples)]
  colnames(rna) <- rna_samples[!duplicated(rna_samples)]
  
  underoverexpressed <- apply(data.matrix(rna), 1, function(x){
    ql <- quantile(as.numeric(unlist(x)), 0.05)
    qu <- quantile(as.numeric(unlist(x)), 0.95)
    x < ql | x > qu
  })
  colnames(underoverexpressed) <- rna_genes
  rownames(underoverexpressed) <- colnames(rna)
  
  data <- list(
    clinical = clinical,
    alteration = data.frame(alt),
    underoverexpressed = data.frame(underoverexpressed)
  )
  saveRDS(data, paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/survival_", C, ".RData"))
}