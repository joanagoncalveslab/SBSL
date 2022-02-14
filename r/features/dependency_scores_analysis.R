library(readr)
library(readxl)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")
source("~/repos/SBSL-modelling-and-analysis/r/utils/global-vars.R")

d2_tissue_types <- c("BREAST", "LUNG", "KIDNEY", "OVARY", "LARGE_INTESTINE", "MULTIPLE_MYELOMA")
avana_tissue_types <- c("BREAST", "LUNG", "KIDNEY", "OVARY", "COLORECTAL", "MULTIPLE_MYELOMA")

all_labels <- labels.load(labels_source)

# load cell line info
sample_info <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/avana/sample_info.csv", 
                        col_types = cols_only(`CCLE Name` = col_guess(), 
                            DepMap_ID = col_guess(), lineage = col_guess()))
# avana dataset
avana <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/avana/Achilles_gene_effect.csv")
avana_genes <- sapply(colnames(avana), function (x) strsplit(x, " \\(")[[1]][1])
avana_samples <- avana$X1
avana$X1 <- NULL
avana <- t(avana)
rownames(avana) <- avana_genes[-1]
colnames(avana) <- avana_samples

# d2 dependency scores
D2 <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/combined_rnai/D2_combined_gene_dep_scores.csv")
d2_genes <- unname(sapply(D2$X1, function (x) strsplit(x, " \\(")[[1]][1]))
d2_samples <- unname(sapply(colnames(D2), function (x) sample_info[sample_info$`CCLE Name` == x, ]$DepMap_ID))
D2$X1 <- NULL
D2 <- as.matrix(D2)
rownames(D2) <- d2_genes
colnames(D2) <- d2_samples[-1]


# ccle mutations
ccle <- read_csv("~/repos/SBSL-modelling-and-analysis/raw_data/depmap/avana/CCLE_mutations.csv", 
                           col_types = cols_only(DepMap_ID = col_guess(), 
                                                 Hugo_Symbol = col_guess(), Tumor_Sample_Barcode = col_guess(), 
                                                 Variant_Classification = col_guess()))
ccle <- ccle[ccle$Variant_Classification != "Silent", ]

calculate_co_dependency <- function(data, geneA, geneB) {
  tryCatch({
    cor.test(as.numeric(data[geneA, ]), as.numeric(data[geneB, ]))
  }
  , error = function(x) c(p.value = 1, estimate = 0)
  , warning = function (x) c(p.value = 1, estimate = 0))
}

calculate_avg_dependency <- function(data, geneA, geneB) {
  tryCatch({
    mean(c(mean(as.numeric(data[geneA, ]), na.rm = TRUE), mean(as.numeric(data[geneB, ]), na.rm = TRUE)))
  }
  , error = function(x) NA
  , warning = function (x) NA)
}

calculate_dependency <- function(data, ccle, geneA, geneB) {
  mutA_samples <- unique(unname(unlist(ccle[which(ccle$Hugo_Symbol == geneA), ][,3])))
  data_mutA <- intersect(mutA_samples, colnames(data))
  out <- c(p.value = 1, statistic = 0)
  
  if (length(data_mutA) >= 3) {
    out <- tryCatch({
      B_scores_mut <- data[geneB, data_mutA]
      B_scores_nonmut <- data[geneB, !(colnames(data) %in% data_mutA)]
      wilcox.test(B_scores_mut, B_scores_nonmut)
    }
    , error = function(x) c(p.value = 1, statistic = 0)
    , warning = function(x) c(p.value = 1, statistic = 0))
  } 
  return(out)
}

for (z in 1:length(cancer_types)) {
  avana_t <- avana_tissue_types[z]
  d2_t <- d2_tissue_types[z]
  C <- cancer_types[z]
  
  # select genes of interest from labels
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  
  # avana samples of interest
  c_samples <- sample_info[sample_info$lineage == tolower(avana_t), ]
  c_avana <- avana[, intersect(avana_samples, unlist(c_samples[,1]))]
  c_d2 <- D2[, intersect(d2_samples, unlist(c_samples[,1]))]
  
  avana_co_dependencies <- numeric(nrow(labels))
  avana_co_dependencies_pvalues <- numeric(nrow(labels))
  avana_avg_dependencies <- numeric(nrow(labels))
  avana_dependencies <- numeric(nrow(labels))
  avana_dependencies_pvalues <- numeric(nrow(labels))
  d2_co_dependencies <- numeric(nrow(labels))
  d2_co_dependencies_pvalues <- numeric(nrow(labels))
  d2_avg_dependencies <- numeric(nrow(labels))
  d2_dependencies <- numeric(nrow(labels))
  d2_dependencies_pvalues <- numeric(nrow(labels))
  for (i in 1:nrow(labels)){
    geneA <- labels[[i, "gene1"]]
    geneB <- labels[[i, "gene2"]]
    
    # get avana dependencies and co-dependencies
    codep <- calculate_co_dependency(c_avana, geneA, geneB)
    avana_co_dependencies[[i]] <- codep["estimate"][[1]]
    avana_co_dependencies_pvalues[[i]] <- codep[["p.value"]]
    avana_avg_dependencies[[i]] <- calculate_avg_dependency(c_avana, geneA, geneB)
    dep <- calculate_dependency(c_avana, ccle, geneA, geneB)
    avana_dependencies[[i]] <- dep["statistic"][[1]]
    avana_dependencies_pvalues[[i]] <- dep[["p.value"]]
    
    # get d2 dependencies and co-dependencies
    codep <- calculate_co_dependency(c_d2, geneA, geneB)
    d2_co_dependencies[[i]] <- codep["estimate"][[1]]
    d2_co_dependencies_pvalues[[i]] <- codep[["p.value"]]
    d2_avg_dependencies[[i]] <- calculate_avg_dependency(c_d2, geneA, geneB)
    dep <- calculate_dependency(c_d2, ccle, geneA, geneB)
    d2_dependencies[[i]] <- dep["statistic"][[1]]
    d2_dependencies_pvalues[[i]] <- dep[["p.value"]]

    if (i%%100 == 0) print(paste("completed", i, "of", nrow(labels)))
  }
  
  results <- data.frame(avana_codep = avana_co_dependencies
                        , avana_codep_pvalue = avana_co_dependencies_pvalues
                        , avana_dep = avana_dependencies
                        , avana_dep_pvalue = avana_dependencies_pvalues
                        , avana_avg = avana_avg_dependencies
                        , d2_codep = d2_co_dependencies
                        , d2_codep_pvalue = d2_co_dependencies_pvalues
                        , d2_dep = d2_dependencies
                        , d2_dep_pvalue = d2_dependencies_pvalues
                        , d2_avg = d2_avg_dependencies)
  
  write.table(results, paste0(output_dir(), C, "_dependency_scores.txt"), sep = "\t", col.names=colnames(results), row.names = FALSE)
}
