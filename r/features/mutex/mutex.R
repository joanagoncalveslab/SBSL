### https://github.com/pathwayanddataanalysis/mutex

library(readr)
library(dplyr)
library(discover)

source("~/repos/msc-thesis-project/r/utils/load-labels.R")
source("~/repos/msc-thesis-project/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types){
  ##########################################################
  ## Mutation Matrix Creation
  ##########################################################
  
  cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", C, "-all_thresholded.by_genes.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", C, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  
  data.mut <- as.matrix(cnv[, -(1:3)])
  rownames(data.mut) <- cnv$`Gene Symbol`
  colnames(data.mut) <- substr(colnames(data.mut), 1, 15)
  
  data.amps <- (data.mut >= 2) * 2
  data.dels <- (data.mut <= -2) * 3
  
  data.mut <- data.amps + data.dels
  
  mutated_genes <- maf[maf$Variant_Classification != "Silent", c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  mutated_genes$Tumor_Sample_Barcode <- substr(mutated_genes$Tumor_Sample_Barcode, 1, 15)
  mutated_genes <- mutated_genes[mutated_genes$Hugo_Symbol %in% cnv$`Gene Symbol`, ]
  mutated_genes <- mutated_genes[mutated_genes$Tumor_Sample_Barcode %in% colnames(data.mut), ]
  mutated_genes <- mutated_genes[!duplicated(mutated_genes), ]
  
  
  for (i in 1:nrow(mutated_genes)) {
    gene <- mutated_genes$Hugo_Symbol[i]
    sample <- mutated_genes$Tumor_Sample_Barcode[i]
    
    if (data.mut[gene, sample] == 0) {
      data.mut[gene, sample] = 1
    } else {
      data.mut[gene, sample] <- data.mut[gene, sample] + 2
    }
    
    if (i%%1000==0) print(paste(i, "of", nrow(mutated_genes)))
  }
  
  ##########################################################
  ## Save file and run mutex java code
  ## Requires parameters file, and for jar file to be precompiled. See git directory for details
  ##########################################################
  
  data.mut <- cbind(Symbol = rownames(data.mut), data.mut)
  dir <- "~/repos/msc-thesis-project/tmp"
  data_file <- paste0(dir, "/dataset.txt")
  write.table(data.mut, data_file, sep = "\t", col.names=colnames(data.mut), row.names=F, quote = FALSE)
  
  
  labels <- all_labels[all_labels$cancer_type == C, ]
  genes <- as.matrix(labels[c("gene1", "gene2")])
  genes <- union(genes[,1], genes[,2])
  genes <- intersect(genes, data.mut[,1])
  genes_file <- paste0(dir, "/genes.txt")
  write.table(genes, genes_file, sep = "\t", row.names=F, col.names = F, quote = FALSE)
  
  system(paste("cp ~/repos/msc-thesis-project/r/features/mutex/parameters.txt", dir))
  
  
  jar_file <- "~/repos/mutex/target/mutex.jar"
  system(paste("java -jar", jar_file, dir))
  
  ##########################################################
  ## Load results
  ##########################################################
  ranked_groups <- read_delim("~/repos/msc-thesis-project/tmp/ranked-groups.txt", 
                              "\t", escape_double = FALSE, col_names = c("score", "gene1", "gene2"), 
                              trim_ws = TRUE, skip = 1)
  
  find_group_score <- function(x) {
    score1 <- ranked_groups[ranked_groups$gene1 == x["gene1"] & ranked_groups$gene2 == x["gene2"], ]$score
    score2 <- ranked_groups[ranked_groups$gene1 == x["gene2"] & ranked_groups$gene2 == x["gene1"], ]$score
    ifelse(length(score1) == 0, score2, score1)
  }
  
  labels$mutex_score <- apply(labels, 1, find_group_score)
  labels$mutex_score[is.na(labels$mutex_score)] <- 2
  
  
  write.table(labels$mutex_score, paste0(output_dir(), C, "_mutex.txt"),
              sep = "\t", col.names=c("mutex"), row.names = FALSE)
  
  print(paste("finished ", C, " analysis"))
}
