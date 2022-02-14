library(discover)
library(survival)
library(survminer)
library(doParallel)
library(foreach)
library(plyr)
library(dplyr)
library(data.table)
library(edgeR)
library(readxl)

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")

# Generate Differential Expression 
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
get_diff_exp <- function (gene1, C) {
  maf_dir <- "~/repos/diff_exp_analysis/maf"
  tumor_feature_counts_file <- "~/repos/diff_exp_analysis/feature_counts/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt"
  manifest_file <- "gdc_manifest.2019-11-28.txt"

  # get cancer types
  GSE62944_06_01_15_TCGA_24_CancerType_Samples <- read.delim(
    "~/repos/diff_exp_analysis/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt", header=FALSE)
  colnames(GSE62944_06_01_15_TCGA_24_CancerType_Samples)[1] <- "sample"
  colnames(GSE62944_06_01_15_TCGA_24_CancerType_Samples)[2] <- "cancer_type"
  samples <- dplyr::filter(GSE62944_06_01_15_TCGA_24_CancerType_Samples, grepl(C, cancer_type))
  # samples <- dplyr::mutate(samples, sample = substr(sample, 1, 16))
  print("loading sample cancer types")

  ## load maf file metadata
  manifest <- read.delim(paste(maf_dir, manifest_file, sep = "/"), header=TRUE)
  print("loaded manifest")
  ## get BRCA samples
  file <- dplyr::filter(manifest, grepl(C, filename))
  maf <- read.delim(paste(maf_dir, file$id, file$filename, sep = "/"), header=TRUE, skip = 5)
  maf <- select(maf, Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification)

  # get mutations only for gene of interest
  maf <- dplyr::filter(maf, Hugo_Symbol == gene1)
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
    et <- exactTest(y)

    # topTags??
    return(et$table)
  }

  return(NULL)
}

# Generate Dependency Scores using any data source
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @param data Matrix in the form of genes * samples. 
# @param mutations MAF matrix with data on mutations for the relevant tissue
# @return Dataframe with all calculated predictors
# 
generate_dependency_scores <- function(gene1, genes, cancer, data, mutations, samples) {
  avgs <- rep(NA, length(genes))
  cor_pvalue <- numeric(length(genes)) + 1
  cor_stat <- numeric(length(genes))
  dep_pvalue <- numeric(length(genes)) + 1
  dep_stat <- numeric(length(genes))
  names(avgs) <- names(cor_pvalue) <- names(cor_stat) <- names(dep_pvalue) <- names(dep_stat) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
    if (!gene2 %in% rownames(data) | !gene1 %in% rownames(data)) next
    
    # Remove NAs
    valuesA <- data[gene1, samples]
    valuesB <- data[gene2, samples]
    nasA <- which(is.na(valuesA))
    nasB <- which(is.na(valuesB))
    nas <- union(nasA, nasB)
    if (length(nas) > 0) {
      valuesA <- valuesA[-nas]
      valuesB <- valuesB[-nas]      
    }
    # Calculate mean
    avgs[gene2] <- mean(c(mean(valuesA), mean(valuesB)))
    # Calculate correlation
    try({
      cor_v <- cor.test(valuesA, valuesB)
      cor_pvalue[gene2] <- cor_v$p.value
      cor_stat[gene2] <- cor_v$estimate
    })

    # Calculate diff dependency for mutation in A
    mutA <- mutations$DepMap_ID[mutations$Hugo_Symbol == gene1]
    scoresB_mutA <- data[gene2, mutA]
    scoresB_nonmutA <- data[gene2, !(colnames(data) %in% mutA)]
    try({
      dep_v <- wilcox.test(scoresB_mutA, scoresB_nonmutA)
      dep_pvalue[gene2] <- dep_v$p.value
      dep_stat[gene2] <- dep_v$statistic
    })
    
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  
  data.frame(
    gene1,     
    gene2 = genes,
    avgs,
    cor_pvalue,
    cor_stat,
    dep_pvalue,
    dep_stat)
}

# Generate CRISPR Dependency Scores
# Data pulled from CRISPR Q3 2019
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_CRISPR_depenency_scores("POLQ", genes, "BRCA")
generate_CRISPR_depenency_scores <- function(gene1, genes, cancer) {
  print(paste0("generating CRISPR dependency scores for ", gene1, " in ", cancer))
  tissues <- c(
    BRCA = "breast",
    LUAD = "lung",
    KIRC = "kidney",
    OV = "ovary",
    COAD = "colorectal",
    LAML = "multiple_myeloma"
  )
  crispr <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/avana.RData")
  sample_info <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/sample_info.RData")
  samples <- colnames(crispr)
  samples <- intersect(samples, sample_info$DepMap_ID[grepl(tissues[cancer], sample_info$lineage)])  
  mutations <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/CCLE_mutations.RData")
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  mutations <- mutations[mutations$DepMap_ID %in% samples, ]
  
  df <- generate_dependency_scores(gene1, genes, cancer, crispr, mutations, samples)
  colnames(df)[3:7] <- c("CRISPR_avg", "CRISPR_cor_pvalue", "CRISPR_cor_stat", "CRISPR_dep_pvalue", "CRISPR_dep_stat")
  df
}

# Generate RNAi Dependency Scores
# Data pulled from RNAi Q3 2019
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_RNAi_depenency_scores("POLQ", genes, "BRCA")
generate_RNAi_depenency_scores <- function(gene1, genes, cancer) {
  print(paste0("generating RNAi dependency scores for ", gene1, " in ", cancer))
  tissues <- c(
    BRCA = "breast",
    LUAD = "lung",
    KIRC = "kidney",
    OV = "ovary",
    COAD = "colorectal",
    LAML = "multiple_myeloma"
  )
  RNAi <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/RNAi_combined.RData")
  sample_info <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/sample_info.RData")
  samples <- colnames(RNAi)
  samples <- intersect(samples, sample_info$DepMap_ID[grepl(tissues[cancer], sample_info$lineage)])  
  mutations <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/CCLE_mutations.RData")
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  mutations <- mutations[mutations$DepMap_ID %in% samples, ]
  
  df <- generate_dependency_scores(gene1, genes, cancer, RNAi, mutations, samples)
  colnames(df)[3:7] <- c("RNAi_avg", "RNAi_cor_pvalue", "RNAi_cor_stat", "RNAi_dep_pvalue", "RNAi_dep_stat")
  df
}


# Generate Gene Expression Correlation in TCGA Patient Tumours
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_gene_expression_correlations("POLQ", genes, "BRCA")
generate_gene_expression_correlations <- function(gene1, genes, cancer) {
  print(paste0("generating gene expression correlations for ", gene1, " in ", cancer))
  mRNA <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/tcga_expression_", cancer, ".RData"))
  mRNA <- as.data.frame(mRNA)
  normal_samples <- substr(colnames(mRNA), 14, 15) == "11"
  tumour_corr <- numeric(length(genes)) 
  tumour_corr.pvalue <- numeric(length(genes)) + 1
  normal_corr <- numeric(length(genes)) 
  normal_corr.pvalue <- numeric(length(genes)) + 1
  names(tumour_corr) <- names(tumour_corr.pvalue) <- names(normal_corr) <- names(normal_corr.pvalue) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
      tumour_valuesA <- as.numeric(mRNA[gene1, !normal_samples])
      tumour_valuesB <- as.numeric(mRNA[gene2, !normal_samples])
    if(!all(is.na(tumour_valuesA)) & !all(is.na(tumour_valuesB))) {
      t.corr <- cor.test(tumour_valuesA, tumour_valuesB, method = "pearson")
      tumour_corr[gene2] <- t.corr[["estimate"]]
      tumour_corr.pvalue[gene2] <- t.corr[["p.value"]]
    }
    
      norm_valuesA <- as.numeric(mRNA[gene1, normal_samples])
      norm_valuesB <- as.numeric(mRNA[gene2, normal_samples])
    if(!all(is.na(norm_valuesA)) & !all(is.na(norm_valuesB))) {
      n.corr <- cor.test(norm_valuesA, norm_valuesB, method = "pearson")
      normal_corr[gene2] <- n.corr[["estimate"]]
      normal_corr.pvalue[gene2] <- n.corr[["p.value"]]    
    }
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(gene1,     
             gene2 = genes,
             tumour_corr,
             tumour_corr.pvalue,
             normal_corr,
             normal_corr.pvalue)
}

# Generate Gene Expression Correlation in GTEx Tissue Samples
# Using v8 data
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_GTEx_scores("POLQ", genes, "BRCA")
generate_GTEx_scores <- function(gene1, genes, cancer) {
  print(paste0("generating GTEx scores for ", gene1, " in ", cancer))
  mRNA <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/GTEx_", cancer, ".RData"))
  gtex_corr <- numeric(length(genes))
  gtex_corr.pvalue <- numeric(length(genes)) + 1
  names(gtex_corr) <- names(gtex_corr.pvalue) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
    valuesA <- as.numeric(mRNA[gene1, ])
    valuesB <- as.numeric(mRNA[gene2, ])
    try({
      n.corr <- cor.test(valuesA, valuesB, method = "pearson")
      gtex_corr[gene2] <- n.corr[["estimate"]]
      gtex_corr.pvalue[gene2] <- n.corr[["p.value"]]
    })
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(gene1,     
             gene2 = genes,
             gtex_corr,
             gtex_corr.pvalue)
}

# Generate Gene Expression Correlation in GTEx Tissue Samples
# Using v8 data
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_pathway_scores("POLQ", genes)
generate_pathway_scores <- function(gene1, genes) {
  print(paste0("generating pathway coparticipation for ", gene1))
  pathways <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/datasets/pathway_gene_sets.RData")
  pathway_coparticipation <- numeric(length(genes)) + 1
  names(pathway_coparticipation) <- genes
  N <- length(pathways)
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
    gene1_pathway <- sapply(pathways, function(x) gene1 %in% x)
    gene2_pathway <- sapply(pathways, function(x) gene2 %in% x)
    K <- sum(gene1_pathway)
    n <- sum(gene2_pathway)
    k <- sum(gene1_pathway & gene2_pathway)
    pathway_coparticipation[gene2] <- phyper(k - 1, K, N - K, n, lower.tail = FALSE, log.p = FALSE)
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(
    gene1,     
    gene2 = genes,
    pathway_coparticipation
  )
}

# Generate Mutex Alteration Scores
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_alteration_mutex_scores("POLQ", genes, "BRCA")
generate_alteration_mutex_scores <- function(gene1, genes, cancer) {
  print(paste0("generating alteration mutex scores for ", gene1, " in ", cancer))
  a <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/mutex_alt_", cancer, ".RData"))
  mutex_alt <- numeric(length(genes)) + 1
  names(mutex_alt) <- genes
  for (gene2 in genes) {
    if (!gene1 %in% rownames(a) | !gene2 %in% rownames(a)) next
    gene1_alt <- a[gene1, ]
    gene2_alt <- a[gene2, ] 
    K <- sum(gene1_alt)
    n <- sum(gene2_alt)
    k <- sum(gene1_alt & gene2_alt)
    N <- ncol(a)
    # hypergeometic test
    mutex_alt[gene2] <- 1 - phyper(k - 1, K, N - K, n, lower.tail = FALSE, log.p = FALSE)
  }
  data.frame(
    gene1,     
    gene2 = genes,
    mutex_alt
  )
}

# Generate Mutex Alteration Scores
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
# Discover Code implemented here: https://github.com/NKI-CCB/DISCOVER
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_alteration_mutex_scores("POLQ", genes, "BRCA")
generate_discover_mutex_scores <- function(gene1, genes, cancer) {
  events <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/discover_event_matrix_", cancer, ".RData")) 
  discover_mutex <- numeric(length(genes)) + 1
  names(discover_mutex) <- genes
  for (gene2 in genes){
    if (!gene1 %in% rownames(events) | !gene2 %in% rownames(events)) next
    result.mutex <- pairwise.discover.test(events[c(gene1, gene2), ])
    discover_mutex[gene2] <- result.mutex$p.values[gene2, gene1]
  }
  data.frame(
    gene1,     
    gene2 = genes,
    discover_mutex
  )
}

# Generate DiscoverSL Mutex Scores
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_alteration_mutex_scores("POLQ", genes, "BRCA")
generate_discoverSL_mutex_scores <- function(gene1, genes, cancer) {
  print(paste0("generating discoverSL mutex scores for ", gene1, " in ", cancer))
  discoverSL.mutex <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/discoversl_mutex_", cancer, ".RData")) 
  scores <- c("amps", "dels", "muts")
  discoversl_mutex <- list()
  for (s in scores) {
    results <- numeric(length(genes)) + 1
    names(results) <- genes
    for (gene2 in genes) {
      valuesA <- tryCatch({discoverSL.mutex[[s]][gene1, ]}, error = function(e) 0)
      valuesB <- tryCatch({discoverSL.mutex[[s]][gene2, ]}, error = function(e) 0)
      K <- sum(valuesA)
      n <- sum(valuesB)
      k <- sum(valuesA & valuesB)
      N <- length(valuesA)
      results[gene2] <- 1 - phyper(k - 1, K, N-K, n, lower.tail = FALSE)
    }
    discoversl_mutex[[s]] <- results
  }
  results <- numeric(length(genes))
  names(results) <- genes
  for (gene2 in genes) {
    results[gene2] <- tryCatch({
      metap::sumlog(c(discoversl_mutex$amps[gene2], discoversl_mutex$dels[gene2], discoversl_mutex$muts[gene2]))$p
    }
    , error = function (e) 1
    , warning = function (e) 1)
  }
  discoversl_mutex$all <- results
  data.frame(
    gene1,     
    gene2 = genes,
    discoversl_mutex_amp = discoversl_mutex$amps,
    discoversl_mutex_del = discoversl_mutex$dels,
    discoversl_mutex_mut = discoversl_mutex$muts,
    discoversl_mutex = discoversl_mutex$all
  )
}

# Generate Survival Mutation Analysis
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_survival_mutation_analysis("POLQ", genes, "BRCA")
generate_survival_mutation_analysis <- function(gene1, genes, cancer) {
  print(paste0("generating survival mutation analysis for ", gene1, " in ", cancer))
  clin <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/clinical_mutation_", cancer, ".RData"))
  # Ugly hack required to get around an issue in survminer and surv_pvalue
  # See https://stackoverflow.com/questions/58991821/error-in-evalinp-data-env-object-database-not-found
  y1_var <- "days_to_death"
  y2_var <- "vital_status"
  x_var <- "test_pair"
  FORMULA <- as.formula(paste("Surv(as.numeric(", y1_var, "), as.numeric(", y2_var, ")) ~ ", x_var))
  mut_logrank_pvals <- numeric(length(genes)) + 1
  names(mut_logrank_pvals) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
    mut_logrank_pvals[gene2] <- tryCatch({
      clin$test_pair <- clin[gene1] & clin[gene2]
      sfit <- survfit(FORMULA, data = clin)
      sfit$call$formula = FORMULA
      p <- as.numeric(surv_pvalue(sfit, data=clin)$pval)
      ifelse(is.na(p), 1, p)
    }, error = function(e) 1)
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(
    gene1,
    gene2 = genes,
    mut_logrank_pvals
  )
}

# Generate Survival Underexpressed mRNA Analysis
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1)
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_survival_underexpressed_mRNA_analysis("POLQ", genes, "BRCA")
generate_survival_underexpressed_mRNA_analysis <- function(gene1, genes, cancer) {
  print(paste0("generating survival mRNA underexpression analysis for ", gene1, " in ", cancer))
  clin <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/clinical_underexpressed_mRNA_", cancer, ".RData"))
  y1_var <- "days_to_death"
  y2_var <- "vital_status"
  x_var <- "test_pair"
  FORMULA <- as.formula(paste("Surv(as.numeric(", y1_var, "), as.numeric(", y2_var, ")) ~ ", x_var))
  mRNA_logrank_pvals <- numeric(length(genes)) + 1
  names(mRNA_logrank_pvals) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  for (gene2 in genes) {
    mRNA_logrank_pvals[gene2] <- tryCatch({
      clin$test_pair <- clin[gene1] & clin[gene2]
      sfit <- survfit(FORMULA, data = clin)
      sfit$call$formula = FORMULA
      p <- as.numeric(surv_pvalue(sfit, data=clin)$pval)
      ifelse(is.na(p), 1, p)
    }, error = function(e) 1)
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(
    gene1,
    gene2 = genes,
    mRNA_logrank_pvals
  )
}

generate_survival_analysis <- function(gene1, genes, cancer) {
  print(paste0("generating survival analysis for ", gene1, " in ", cancer))
  survival_data <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/survival_", cancer, ".RData"))
  y1_var <- "days_to_death"
  y2_var <- "vital_status"
  x_var <- "test_pair + gender + years_to_birth + race"
  FORMULA <- as.formula(paste("Surv(as.numeric(", y1_var, "), as.numeric(", y2_var, ")) ~ ", x_var))
  logrank_pvals <- numeric(length(genes)) + 1
  names(logrank_pvals) <- genes
  i <- 0
  pb <- txtProgressBar(min = i, max = length(genes), style = 3)
  
  clin <- survival_data$clinical
  clin$years_to_birth <- as.numeric(as.character(clin$years_to_birth))
  for (gene2 in genes) {
    logrank_pvals[gene2] <- tryCatch({
      gene1_inactive <- survival_data$alteration[row.names(clin), gene1] | survival_data$underoverexpressed[row.names(clin), gene1]
      gene1_inactive[is.na(gene1_inactive)] <- FALSE
      gene2_inactive <- survival_data$alteration[row.names(clin), gene2] | survival_data$underoverexpressed[row.names(clin), gene2]
      gene2_inactive[is.na(gene2_inactive)] <- FALSE
      clin$test_pair <- gene1_inactive & gene2_inactive
      res.cox <- coxph(FORMULA, data = clin)
      p <- summary(res.cox)$coefficients[1, 5]
      ifelse(is.na(p), 1, p)
    }, error = function(e) 1)
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  data.frame(
    gene1,
    gene2 = genes,
    logrank_pvals
  )
}

# Generate mutual exclusivity scores based on the MUTEX algorithm
# Data pulled from RTCGAToolbox::getFirehoseRunningDates(1) and mutex scores generated using code from
# https://github.com/pathwayanddataanalysis/mutex
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_MUTEX_scores("POLQ", genes, "BRCA")
generate_MUTEX_scores <- function(gene1, genes, cancer) {
  MUTEX <- numeric(length(genes)) + 2
  names(MUTEX) <- genes
  try({
    ranked_groups <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/datasets/MUTEX_", cancer, ".RData"))
    for (gene2 in genes) {
      if (!gene1 %in% rownames(ranked_groups) | !gene2 %in% rownames(ranked_groups)) next
      score1 <- ranked_groups[ranked_groups$gene1 == gene1 & ranked_groups$gene2 == gene2, ]$score
      score2 <- ranked_groups[ranked_groups$gene1 == gene2 & ranked_groups$gene2 == skyoe, ]$score
      MUTEX[gene2] <- ifelse(length(score1) == 0, score2, score1)
    }
  })
  MUTEX[is.na(MUTEX)] <- 2
  data.frame(
    gene1,
    gene2 = genes,
    MUTEX
  )
}

# Generate differential expression scores
# Data pulled from GEO
#
# @param gene1 The query gene
# @param cancer The cancer of interest
# @return Dataframe with all calculated predictors
# 
# @example generate_diff_expression_scores("POLQ", genes, "BRCA")
generate_diff_expression_scores <- function(gene1, genes, cancer) {
  diff_exp <- get_diff_exp(gene1, cancer)
  diff_exp_logFC <- numeric(length(genes))
  diff_exp_pvalue <- numeric(length(genes)) + 1
  names(diff_exp_logFC) <- names(diff_exp_pvalue) <- genes
  for (gene2 in genes) {    
    pvalue1 <- diff_exp[gene2, ]$PValue
    diff_exp_pvalue[gene2] <- ifelse(length(pvalue1), pvalue1, 1) 
    
    logFC1 <- diff_exp[gene2, ]$logFC
    diff_exp_logFC[gene2] <- ifelse(length(pvalue1), logFC1, 0)
  }
  data.frame(
    gene1,
    gene2 = genes,
    diff_exp_logFC,
    diff_exp_pvalue
  )
}


generate_all_scores_for_gene <- function(gene1, genes, cancer, SL) {
  # details
  df <- data.frame(
   gene1,
   gene2 = genes,
   cancer_type = cancer,
   SL
  )
  # survival
  df1 <- generate_survival_underexpressed_mRNA_analysis(gene1, genes, cancer)
  df2 <- generate_survival_mutation_analysis(gene1, genes, cancer)
  # # mutex
  df3 <- generate_discoverSL_mutex_scores(gene1, genes, cancer)
  df4 <- generate_discover_mutex_scores(gene1, genes, cancer)
  df5 <- generate_MUTEX_scores(gene1, genes, cancer)
  df6 <- generate_alteration_mutex_scores(gene1, genes, cancer)
  # # dependencies
  df7 <- generate_RNAi_depenency_scores(gene1, genes, cancer)
  df8 <- generate_CRISPR_depenency_scores(gene1, genes, cancer)
  # # expression correlation
  df9 <- generate_GTEx_scores(gene1, genes, cancer)
  df10 <- generate_gene_expression_correlations(gene1, genes, cancer)
  # # other
  df11 <- generate_diff_expression_scores(gene1, genes, cancer)
  df12 <- generate_pathway_scores(gene1, genes)
  df13 <- generate_survival_analysis(gene1, genes, cancer)
  
  combined_df <- Reduce(
    function(x, y, ...) merge(x, y, by=c("gene1", "gene2")), 
    list(df, df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13)
  )
  combined_df
}

generate_all_scores_for_dataset <- function(dataset, cancer_list) {
  num_cores <- parallel::detectCores()
  registerDoParallel(cores=num_cores)
  pairs <- labels.load(dataset)
  
  pairs <- pairs[pairs$cancer_type %in% cancer_list, ]
  
  if ("LAML" == cancer_list) {
    pos_samples <- pairs[pairs$SL == 1, ]
    neg_samples <- pairs[pairs$SL == 0, ]
    set.seed(123)
    n_pos_samples <- nrow(pos_samples)
    n_neg_samples <- nrow(neg_samples)
    n_samples <- min(n_pos_samples, n_neg_samples)
    pos_samples <- dplyr::sample_n(pos_samples, n_samples)
    neg_samples <- dplyr::sample_n(neg_samples, n_samples)
    pairs <- dplyr::bind_rows(list(pos_samples, neg_samples))
  }
  # pairs <- readRDS("~/repos/SBSL-modelling-and-analysis/r/data/isle_new_LAML.Rdata")[1:4]
  gene1s <- unique(unlist(pairs$gene1))
  
  r <- foreach(g = 1:length(gene1s)) %dopar% {
    l <- list()
    gene1 <- gene1s[g]
    d <- pairs[pairs$gene1 == gene1, ]
    cancer_types <- unique(d$cancer_type)
    for (cancer in cancer_types) {
      d_ <- d[d$cancer_type == cancer, ]
      non_dups <- which(!duplicated(d_$gene2))
      genes <- d_$gene2[non_dups]
      SL <- d_$SL[non_dups]
      l[[cancer]] <- generate_all_scores_for_gene(gene1, genes, cancer, SL)
    }
    print(paste(g, "of", length(gene1s), "complete"))
    dplyr::bind_rows(l)
  }
  dplyr::bind_rows(r)
}

# Example for single gene
# gene1 <- "POLQ"
# genes <- c("BRCA1", "BRCA2", "PARP1")
# cancer <- "BRCA"
# cancer_types <- c("BRCA", "LUAD", "KIRC", "OV")
# generate_all_scores_for_gene(gene1, genes, cancer, c(1,1,1))

# Example for dataset
# dataset <- "isle"
# cancer_list <- c("LAML")
# laml <- generate_all_scores_for_dataset(dataset, cancer_list)
# saveRDS(laml, "~/repos/SBSL-modelling-and-analysis/r/data/isle_new_LAML_exp.Rdata")
