library(plyr)
library(dplyr)
library(data.table)
library(edgeR)
library(readxl)

args = commandArgs(trailingOnly=TRUE)

start_index = as.numeric(args[1])
end_index = start_index + as.numeric(args[2])
C <- args[3]

maf_dir <- "./maf"
tumor_feature_counts_file <- "./feature_counts/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt"
manifest_file <- "gdc_manifest.2019-11-28.txt"

print(C)

# get cancer types
GSE62944_06_01_15_TCGA_24_CancerType_Samples <- read.delim(
  "./GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt", header=FALSE)
colnames(GSE62944_06_01_15_TCGA_24_CancerType_Samples)[1] <- "sample"
colnames(GSE62944_06_01_15_TCGA_24_CancerType_Samples)[2] <- "cancer_type"
samples <- dplyr::filter(GSE62944_06_01_15_TCGA_24_CancerType_Samples, grepl(C, cancer_type))
# samples <- dplyr::mutate(samples, sample = substr(sample, 1, 16))
print("loading sample cancer types")

## Load processed gene list
labels <- read_excel("./ISLE.xlsx", skip = 2)
labels <- dplyr::filter(labels, labels$`cancer type tested` == C)
genes <- as.matrix(labels[c("gene1", "gene2")])
genes <- union(genes[,1], genes[,2])

print("loaded genes")
## load maf file metadata
manifest <- read.delim(paste(maf_dir, manifest_file, sep = "/"), header=TRUE)
print("loaded manifest")
## get BRCA samples
file <- dplyr::filter(manifest, grepl(C, filename))
maf <- read.delim(paste(maf_dir, file$id, file$filename, sep = "/"), header=TRUE, skip = 5)
maf <- select(maf, Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification)

# get mutations only for gene of interest
maf <- dplyr::filter(maf, Hugo_Symbol %in% genes)
maf <- dplyr::filter(maf, !grepl("Silent", Variant_Classification))
maf <- dplyr::mutate(maf, Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 16))
print("loaded mutations")

## Load raw count tables
tumor_feature_counts <- read.delim(tumor_feature_counts_file, sep = "\t", header = TRUE)
row.names(tumor_feature_counts) <- tumor_feature_counts$X
tumor_feature_counts$X <- NULL
colnames(tumor_feature_counts) <- gsub("\\.", "-", colnames(tumor_feature_counts))
print("loaded feature counts")
print(dim(tumor_feature_counts))

# select brca only samples
tumor_feature_counts <- select(tumor_feature_counts, samples$sample)
print("subsetting to select only tissue type samples")
print(dim(tumor_feature_counts))

## diffExpression
if (end_index > length(genes)) {
  end_index = length(genes)
}

for (gene1 in genes[start_index:end_index]) {
  print(gene1)
  # select driver gene
  gene1 <- as.character(gene1)
  
  # select samples with mutation in driver gene
  samples_with_gene1_mutation <- dplyr::filter(maf, gene1 == Hugo_Symbol)

  # mark as TRUE any sample which has a mutation in the driver gene, which will create two groups
  groups <- sapply(colnames(tumor_feature_counts), function(x){ substr(x, 1, 16) %in% samples_with_gene1_mutation$Tumor_Sample_Barcode})
  print(sum(groups))

  if (sum(groups) >= 1) {
      # estimate disparity
    y <- DGEList(counts = tumor_feature_counts, group = groups)

    print("normalising")
    # trimmed mean of m values normalisation
    y <- calcNormFactors(y)

    # create design matrix
    design <- model.matrix(~groups)

    print("estimating disparity")
    # estimate disparity
    y <- estimateDisp(y,design)

    print("performing tests")
    # exact test
    et <- exactTest(y[intersect(rownames(y), genes), ])

    # topTags??
    write.table(et$table, file = paste0("./", C, "_results/", gene1, "__", sum(groups)), sep = "\t", col.names=NA, row.names = T)
    print(gene1)
  }
}

print("finished")
