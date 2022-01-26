rm(list=ls())
source("~/repos/msc-thesis-project/r/utils/load-raw.R")
source("~/repos/msc-thesis-project/r/utils/load-labels.R")
cancer <- "LAML"
###
# TCGA Data
###
tcga_expression <- raw.get_tcga_expression(cancer)
saveRDS(tcga_expression, paste0("~/repos/msc-thesis-project/r/data/datasets/tcga_expression_", cancer, ".RData"))
rm(list=ls())
###
# GTEx Data
###
cancer_genes <- labels.load("isle")
cancer_genes <- cancer_genes[cancer_genes$cancer_type == cancer, ]
cancer_genes <- union(cancer_genes$gene1, cancer_genes$gene2)
GTEx_expression <- raw.get_GTEx_expression(cancer, cancer_genes)
saveRDS(GTEx_expression, paste0("~/repos/msc-thesis-project/r/data/datasets/GTEx_", cancer, ".RData"))
rm(list=ls())
###
# Mutation Alteration
###
cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
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
saveRDS(data.mut, paste0("~/repos/msc-thesis-project/r/data/datasets/mutex_alt_", cancer, ".RData"))
rm(list=ls())
###
# Discover Event Matrix
###
library(discover)
cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
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
events <- discover.matrix(data.mut)
saveRDS(events, paste0("~/repos/msc-thesis-project/r/data/datasets/discover_event_matrix_", cancer, ".RData")) 
rm(list=ls())
###
# DiscoverSL
### 
cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
print("cnv and genes loaded")
row.names(cnv) <- cnv$`Gene Symbol`
symbols <- cnv$`Gene Symbol`
cnv <- cnv[,4:dim(cnv)[2]]
row.names(cnv) <- symbols
x1 <- as.matrix(cnv)
discoversl_mut <- list()
discoversl_mut$amps <- x1 >= 2
discoversl_mut$dels <- x1 <= -2
maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
maf <- dplyr::select(maf, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
maf <- dplyr::filter(maf, Variant_Classification != "Silent")
maf <- dplyr::select(maf, Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)
x <- as.matrix(table(maf[c("Hugo_Symbol", "Tumor_Sample_Barcode")])) > 0
discoversl_mut$muts <- x
saveRDS(discoversl_mut, paste0("~/repos/msc-thesis-project/r/data/datasets/discoversl_mutex_", cancer, ".RData")) 
rm(list=ls())
###
# MUTEX
###
### https://github.com/pathwayanddataanalysis/mutex
library(readr)
library(dplyr)
cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
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
ranked_groups <- read_delim("~/repos/msc-thesis-project/tmp/ranked-groups.txt", 
                              "\t", escape_double = FALSE, col_names = c("score", "gene1", "gene2"), 
                              trim_ws = TRUE, skip = 1)
saveRDS(ranked_groups, paste0("~/repos/msc-thesis-project/r/data/datasets/MUTEX_", cancer, ".RData"))
rm(list=ls())
###
# TCGA Survival
###
clinical <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Clinical.txt"), 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
cols <- unname(unlist(clinical[, 1]))
clinical[, 1] <- NULL
clinical <- as.data.frame(t(clinical))
colnames(clinical) <- cols

cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
cnv_genes <- unlist(cnv[, 1])
cnv <- cnv[, -c(1, 2, 3)]
cnv_samples <- lapply(colnames(cnv), function (x) tolower(substr(x, 1, 12)))
cnv <- as.data.frame(t(cnv))
rownames(cnv) <- cnv_samples
colnames(cnv) <- cnv_genes

maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
maf <- dplyr::filter(maf, Variant_Classification != "Silent")
maf <- dplyr::select(maf, Hugo_Symbol, Tumor_Sample_Barcode)
maf$Tumor_Sample_Barcode <- unlist(lapply(maf$Tumor_Sample_Barcode, function (x) tolower(substr(x, 1, 12))))

alt <- (cnv <= -2) | (cnv >= 2)
for (i in 1:nrow(maf)) {
  gene <- maf$Hugo_Symbol[i]
  sample <- maf$Tumor_Sample_Barcode[i]
  tryCatch({
    alt[sample, gene] <- alt[sample, gene] | TRUE
  }, error = function(x){})
  if (i%%100 == 0) print(paste(i, "of", nrow(maf)))
}
rna <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-RNAseq2GeneNorm.txt"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
rna <- rna[-(1:30), ]
rna_samples <- unname(sapply(colnames(rna)[-1], function (x) tolower(substr(x, 1, 12))))
rna_genes <- unname(sapply(unlist(rna[, 1]), function (x) strsplit(x, "\\|")[[1]][1]))
rna <- rna[, -1]
rna <- rna[, !duplicated(rna_samples)]
colnames(rna) <- rna_samples[!duplicated(rna_samples)]

underexpressed <- as.data.frame(apply(rna, 1, function(x) x < quantile(as.numeric(unlist(x)), 0.05)))
colnames(underexpressed) <- rna_genes
rownames(underexpressed) <- colnames(rna)
alt_me <- merge(clinical, alt, by=0)
under_me <- merge(clinical, underexpressed, by=0)
saveRDS(alt_me, paste0("~/repos/msc-thesis-project/r/data/datasets/clinical_mutation_", cancer, ".RData"))
saveRDS(under_me, paste0("~/repos/msc-thesis-project/r/data/datasets/clinical_underexpressed_mRNA_", cancer, ".RData"))
rm(list=ls())