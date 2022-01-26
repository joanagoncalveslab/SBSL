# Following the methods outlined here
# https://www.sciencedirect.com/science/article/pii/S0092867414009775?via%3Dihub

library(readr)
library(dplyr)

source("~/repos/msc-thesis-project/r/utils/load-raw.R")

genomic_survival_of_the_fittest <- function(labels, cna, deleted, mRNA, mut, underexpressed) {
  gsof <- numeric(nrow(labels))
  for (i in  1:nrow(labels)) {
    out_gsof <- numeric(2)
    for (j in 1:2) {
      geneA <- ifelse(j == 1, "gene1", "gene2")
      geneB <- ifelse(geneA == "gene1", "gene2", "gene1")
      geneA <- labels[[i, geneA]]
      geneB <- labels[[i, geneB]]
      inactive_samples <- tryCatch({
        geneA_deleted_samples <- colnames(cna)[unname(deleted[geneA, ])]
        geneA_underexpressed_samples <- colnames(mRNA)[unname(unlist(underexpressed[geneA, ]))]
        geneA_mutated_samples <- mut[mut["Hugo_Symbol"] == geneA, ]["Tumor_Sample_Barcode"]
        unlist(union(union(geneA_mutated_samples, geneA_deleted_samples), geneA_underexpressed_samples))
      }
      ,error = function(e) {
        c()
      })

      out_gsof[j] <- tryCatch({
        x <- as.numeric(cna[geneB, intersect(colnames(cna),inactive_samples)])
        y <- as.numeric(cna[geneB, !colnames(cna) %in% inactive_samples])
        wilcox.test(x, y, alternative = "less")$p.value
      }
      , error = function (e) {
        1
      })
    }
    gsof[i] <- min(out_gsof)
    if(i%%100 == 0) print(paste("gsof", i, "of", nrow(labels)))
  }
  gsof
} 


shRNA_functional_based_examination <- function(labels, cna, deleted, mRNA, mut, underexpressed, shRNA) {
  shRNAbfes <- numeric(nrow(labels))
  for (i in  1:nrow(labels)) {
    out_shRNAbfe <- numeric(2)
    for (j in 1:2) {
      geneA <- ifelse(j == 1, "gene1", "gene2")
      geneB <- ifelse(geneA == "gene1", "gene2", "gene1")
      geneA <- labels[[i, geneA]]
      geneB <- labels[[i, geneB]]
      inactive_samples <- tryCatch({
        geneA_deleted_samples <- colnames(cna)[unname(deleted[geneA, ])]
        geneA_underexpressed_samples <- colnames(mRNA)[unname(unlist(underexpressed[geneA, ]))]
        geneA_mutated_samples <- mut[mut["Hugo_Symbol"] == geneA, "Tumor_Sample_Barcode"]
        unlist(union(union(geneA_mutated_samples, geneA_deleted_samples), geneA_underexpressed_samples))
      }
      ,error = function(e) {
        c()
      })

      out_shRNAbfe[j] <- tryCatch({
        x <- as.numeric(shRNA[geneB, intersect(colnames(shRNA),inactive_samples)])
        y <- as.numeric(shRNA[geneB, !colnames(shRNA) %in% inactive_samples])
        wilcox.test(x, y, alternative = "less")$p.value
      }
      , error = function (e) 1)
    }
    shRNAbfes[i] <- min(out_shRNAbfe)
    if(i%%100 == 0) print(paste("shRNAbfes", i, "of", nrow(labels)))
  }
  shRNAbfes
}

expression_correlation <- function(labels, mRNA) {
  mRNA_corr <- numeric(nrow(labels))
  mRNA_corr.pvalue <- numeric(nrow(labels))
  
  for (i in  1:nrow(labels)) {
    ccle.corr <- tryCatch({
      gene1_exp_tumour <- as.numeric(mRNA[labels[[i, "gene1"]], ])
      gene2_exp_tumour <- as.numeric(mRNA[labels[[i, "gene2"]], ])
      cor.test(gene1_exp_tumour, gene2_exp_tumour, method = "spearman")
    }, error = function(e) {
      c(estimate = 0, p.value = 1)
    }, warning = function(e) {
      c(estimate = 0, p.value = 1)
    })
    
    mRNA_corr[i] <- ccle.corr[["estimate"]]
    mRNA_corr.pvalue[i] <- ccle.corr[["p.value"]]
  }
  return(list(coef = mRNA_corr, pvalue = mRNA_corr.pvalue))
}


DAISY <- function(gene_pairs){
  cancers <- unique(gene_pairs$cancer_type)
  results <- list()
  
  for (c in 1:length(cancers)) {
    C <- cancers[c]
    print(paste("starting ", C, " analysis"))
    labels <- dplyr::filter(gene_pairs, gene_pairs$cancer_type == C)
    genes <- union(labels$gene1, labels$gene2)
    
    #########################################
    ## Set up pipeline
    #########################################
    # thresholds defined in supplementary materials
    ccle_thres_del <- -0.3
    tcga_thres_del <- -0.3
    thres_underexpressed <- .1
    
    ### Load CCLE Data
    shRNA <- raw.get_ccle_rnai_dependency_scores(C, genes)
    cna <- raw.get_ccle_rnai_cna(C, genes)
    deleted <- cna < ccle_thres_del
    mRNA <- raw.get_ccle_rnai_expression(C)
    underexpressed <- apply(mRNA, 2, function(x) x < quantile(as.numeric(unlist(x)), thres_underexpressed, na.rm = TRUE))
    mut <- raw.get_ccle_rnai_mutations(C, genes)
    # Get p-values for CCLE Data
    labels$ccle_sof.pvalue <- genomic_survival_of_the_fittest(labels, cna, deleted, mRNA, mut, underexpressed)
    labels$ccle_shRNA.pvalue <- shRNA_functional_based_examination(labels, cna, deleted, mRNA, mut, underexpressed, shRNA)
    ccle_mRNA_corr <- expression_correlation(labels, mRNA)
    labels$ccle_mRNA.coef <- ccle_mRNA_corr$coef
    labels$ccle_mRNA.pvalue <- ccle_mRNA_corr$pvalue
    
    rm(cna, deleted, mRNA, underexpressed, mut, shRNA)
    
    ### Load TCGA Data
    cna <- raw.get_tcga_cna(C, genes, thresholded = FALSE)
    deleted <- cna <= tcga_thres_del
    mRNA <- raw.get_tcga_expression(C)
    underexpressed <- apply(mRNA, 2, function(x) x < quantile(as.numeric(unlist(x)), thres_underexpressed, na.rm = TRUE))
    mut <- raw.get_tcga_mutations(C, genes)
    # Get p-values for TCGA Data
    labels$tcga_sof.pvalue <- genomic_survival_of_the_fittest(labels, cna, deleted, mRNA, underexpressed, mut)
    tcga_mRNA_corr <- expression_correlation(labels, mRNA)
    labels$tcga_mRNA.coef <- tcga_mRNA_corr$coef
    labels$tcga_mRNA.pvalue <- tcga_mRNA_corr$pvalue
    rm(cna, deleted, mRNA, underexpressed, mut)
    
    print(paste("Spearman Coefficients Calulated"))
    results[[c]] <- labels
    print(paste(C, "finished"))
  }
  
  ##########################################
  ## Combine p-values and do corrections
  ##########################################
  gene_pairs <- bind_rows(results)
  
  ##########################################
  ## Combine results from different datasets
  ##########################################
  gene_pairs$sof.pvalue <- apply(gene_pairs, 1, function(x) tryCatch(
    {metap::sumlog(as.numeric(c(x["ccle_sof.pvalue"], x["tcga_sof.pvalue"])))$p}
    ,error = function(e) NA
  ))
  gene_pairs$sof.pvalue.corrected <- p.adjust(gene_pairs$sof.pvalue, "bonferroni")
  
  
  gene_pairs$mRNA.pvalue <- apply(gene_pairs, 1, function(x) tryCatch(
    {metap::sumlog(as.numeric(c(x["ccle_mRNA.pvalue"], x["tcga_mRNA.pvalue"])))$p}
    ,error = function(e) NA
  ))
  gene_pairs$mRNA.pvalue.corrected <- p.adjust(gene_pairs$mRNA.pvalue, "bonferroni")
  
  ##########################################
  ## Combine into one p-value and perform correction
  ##########################################
  gene_pairs$all.pvalue <- apply(gene_pairs, 1, function(x) tryCatch(
    {metap::sumlog(as.numeric(c(x["sof.pvalue.corrected"], x["mRNA.pvalue.corrected"], x["ccle_shRNA.pvalue"])))$p}
    ,error = function(e) NA
  ))
  gene_pairs$all.pvalue.corrected <- p.adjust(gene_pairs$all.pvalue, "bonferroni")
  
  ###########################################
  ## Mark SL Pairs
  ###########################################
  gene_pairs$daisy_sl_pair <- gene_pairs$all.pvalue.corrected < 0.05
  print("returning results")
  return(gene_pairs)
}

DAISY.predict <- function(test) {
  data_dir <- "~/repos/msc-thesis-project/r/data/"
  daisy_results <- readRDS(paste0(data_dir, "DAISY_results.RData"))
  no_result_default <- 1e-7 # return value if no precalculated results exists
  unname(apply(test, 1, function(x, d) {
    gene1 <- x[[1]]
    gene2 <- x[[2]]
    cancer <- x[[3]]
    matched_rows <- d[d$gene1 == gene1 & d$gene2 == gene2 & d$cancer_type == cancer, ncol(d) - 2]
    if (length(matched_rows) > 0) {
      r <- ifelse(is.na(matched_rows[1]), no_result_default, matched_rows[1])
    }
    ifelse(round(r, 8) == 0, no_result_default, r) # floating point issue fix
  }, d = daisy_results))
}
