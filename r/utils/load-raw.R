library(readr)
library("org.Hs.eg.db")

raw.get_tcga_expression <- function(cancer, genes = NA) {
  rna <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", cancer, "-RNAseq2GeneNorm.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  rna <- rna[2:dim(rna)[1], ]
  colnames(rna)[1] <- "gene_id"
  rna$Hugo_Symbol <- sapply(rna$gene_id, function (x) strsplit(x, "\\|")[[1]][1])
  
  if (typeof(genes) == "character") {
    rna <- dplyr::filter(rna, Hugo_Symbol %in% genes)  
  }
  
  rna <- rna[!duplicated(rna$Hugo_Symbol), ]
  symbols <- rna$Hugo_Symbol
  rna$Hugo_Symbol <- NULL
  rna$gene_id <- NULL
  row.names(rna) <- symbols
  rna
}

raw.get_tcga_mutations <- function(cancer, genes = NA) {
  mutations <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  if (typeof(genes) == "character") {
    mutations <- mutations[mutations$Hugo_Symbol %in% genes, ]
  }
  mutations
}

raw.get_tcga_cna <- function(cancer, genes = NA, thresholded = TRUE) {
  if (thresholded) {
    cnv <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                      "\t", escape_double = FALSE, trim_ws = TRUE)   
  } else {
    cnv <- read_delim(paste0("~/repos/SBSL-modelling-and-analysis/raw_data/firehose/20160128-", cancer, "-all_data_by_genes.txt"), 
                      "\t", escape_double = FALSE, trim_ws = TRUE)    
  }
  cnv_genes <- unlist(cnv$`Gene Symbol`)
  rownames(cnv) <- cnv_genes
  cnv_samples <- colnames(cnv)[4:ncol(cnv)]
  
  if (typeof(genes) == "character") {
    cnv_genes <- intersect(cnv_genes, genes)
  }
  
  result <- cnv[cnv_genes, cnv_samples]
  colnames(result) <- cnv_samples
  rownames(result) <- cnv_genes
  result
}

raw.get_GTEx_expression <- function (cancer, genes) {
  sample_attr <- read_delim("~/repos/SBSL-modelling-and-analysis/raw_data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", 
                            "\t", escape_double = FALSE, col_types = cols_only(SAMPID = col_guess(), 
                                                                               SMTS = col_guess(),
                                                                               SMTSD = col_guess()), trim_ws = TRUE)
  
  cancer2tissue <- c(BRCA = "Breast - Mammary Tissue"
                     , LUAD = "Lung"
                     , KIRC = "Kidney - Cortex"
                     , OV = "Ovary"
                     , COAD = "Colon"
                     , LAML = "Blood")
  
  samples <- sample_attr$SAMPID[grepl(cancer2tissue[cancer], sample_attr$SMTSD)]
  print(paste("loaded", length(samples), cancer, "samples"))
  
  # Get ensembl ids
  ensembl_ids <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENSEMBL")
  genes <- genes[!is.na(ensembl_ids)]
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
  names(genes) <- ensembl_ids
  print("got ensembl ids")
  
  # first fine the correct lines to read
  f <- "~/repos/SBSL-modelling-and-analysis/raw_data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
  conn <- file(f, "r")
  lines <- list()
  names <- read.table(conn, sep = "\t", skip = 2, nrow = 1)
  names <- unname(unlist(names))
  cols_to_keep <- c(1, which(names %in% samples))
  lines_in_file <- 56203 - 3 # we have already read 3 lines by this point
  for (i in 1:lines_in_file) {
    r <- readLines(conn, 1)
    if (substr(r, 1, 15) %in% ensembl_ids) {
      lines[[length(lines) + 1]] <- unlist(strsplit(r, "\t"))[cols_to_keep]
    }
    if (i%%100 == 0) print(paste("completed", i, "of", lines_in_file))
  }
  close(conn)
  
  print("binding rows")
  df <- data.frame(do.call(rbind, lines))
  colnames(df) <- names[cols_to_keep]
  rownames(df) <- unname(sapply(substr(df$Name, 1, 15), function(x) genes[x]))
  df$Name <- NULL
  
  df
}

raw.get_ccle_rnai_dependency_scores <- function(cancer, genes = NULL) {
  sample_info <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/sample_info.csv")
  cancer2tissue <- c(BRCA = "breast"
                     , LUAD = "lung"
                     , KIRC = "kidney"
                     , OV = "ovary"
                     , COAD = "colorectal")
  c_samples <- sample_info[sample_info$disease == cancer2tissue[cancer], ]
  
  D2 <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/D2_combined_gene_dep_scores.csv")
  D2_genes <- sapply(strsplit(D2$X1, " "), function (x) x[[1]])
  D2 <- D2[!duplicated(D2_genes), ]
  D2$X1 <- NULL
  D2 <- D2[, intersect(colnames(D2), unlist(c_samples[,1]))]
  rownames(D2) <- D2_genes[!duplicated(D2_genes)]
  if (!is.null(genes)) {
    common_genes <- intersect(D2_genes, genes)
    D2 <- as.data.frame(D2)[common_genes, intersect(colnames(D2), unlist(c_samples[,1]))]
  }
  as.data.frame(D2)
}

raw.get_ccle_rnai_cna <- function(cancer, genes) {
  sample_info <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/sample_info.csv")
  cancer2tissue <- c(BRCA = "breast"
                     , LUAD = "lung"
                     , KIRC = "kidney"
                     , OV = "ovary"
                     , COAD = "colorectal")
  c_samples <- sample_info[sample_info$disease == cancer2tissue[cancer], ]
  
  cn <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/WES_SNP_CN_data.csv")
  cn_genes <- sapply(strsplit(cn$X1, " "), function (x) x[[1]])
  cn <- cn[!duplicated(cn_genes), ]
  rownames(cn) <- cn_genes[!duplicated(cn_genes)]
  cn$X1 <- NULL
  common_genes <- intersect(cn_genes, genes)
  cn <- cn[common_genes, intersect(colnames(cn), unlist(c_samples[,1]))]
  rownames(cn) <- common_genes
  cn
}

raw.get_ccle_rnai_mutations <- function(cancer, genes) {
  sample_info <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/sample_info.csv")
  cancer2tissue <- c(BRCA = "breast"
                     , LUAD = "lung"
                     , KIRC = "kidney"
                     , OV = "ovary"
                     , COAD = "colorectal")
  c_samples <- sample_info[sample_info$disease == cancer2tissue[cancer], ]
  mutations <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/CCLE_mutation_data.csv")
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  mutations <- mutations[mutations$Hugo_Symbol %in% genes, ]
  mutations <- mutations[mutations$Tumor_Sample_Barcode %in% c_samples$CCLE_ID, ]
  mutations
}

raw.get_ccle_rnai_expression <- function(cancer) {
  sample_info <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/sample_info.csv")
  cancer2tissue <- c(BRCA = "breast"
                     , LUAD = "lung"
                     , KIRC = "kidney"
                     , OV = "ovary"
                     , COAD = "large_intestine")
  c_samples <- sample_info[sample_info$disease == cancer2tissue[cancer], ]
  
  ex <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/RNAseq_lRPKM_data.csv")
  ex_genes <- sapply(strsplit(ex$X1, " "), function (x) x[[1]])
  ex <- ex[!duplicated(ex_genes), ]
  rownames(ex) <- ex_genes[!duplicated(ex_genes)]
  ex$X1 <- NULL
  ex[, intersect(colnames(ex), unlist(c_samples[,1]))]
  rownames(ex) <- ex_genes[!duplicated(ex_genes)]
  ex
}
