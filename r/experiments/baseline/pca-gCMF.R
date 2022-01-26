##################################################
# pca-gCMF_data_generation
#
# Generate data files for pca-gCMF method
##################################################
library("Hmisc")
library(readr)
source("~/repos/msc-thesis-project/r/utils/load-raw.R")
source("~/repos/msc-thesis-project/r/utils/train-model.R")
source("~/repos/msc-thesis-project/r/experiments/baseline/source/pca_gCMF.R")

pca_gCMF.generate_feature_matrices <- function(cancer) {
  d <- train.get_dataset("combined")
  all_genes <- union(d$gene1, d$gene2)
  output <- paste0("~/repos/msc-thesis-project/r/data/pca-gCMF/", cancer, "/")

  ##################################################
  ## RNA expression
  ##################################################
  # get RNA
  exp <- raw.get_tcga_expression(cancer, all_genes)
  genes <- row.names(exp)
  normal_samples <- substr(colnames(exp), 14, 15) == "11"
  exp <- exp[, !normal_samples]
  rownames(exp) <- genes
  # create empty matrix
  data_RNA_expression_alltraintest <- data.frame(matrix(NA, nrow = length(all_genes), ncol = ncol(exp)), row.names = all_genes)
  for (g in genes){
    data_RNA_expression_alltraintest[g, ] <- exp[g, ]
  }

  ##################################################
  ## Pairwise Co-expression
  ##################################################
  coexp <- rcorr(t(data.matrix(exp)), type = "spearman")$P
  coexp <- matrix(p.adjust(coexp, method = "bonferroni"), nrow(coexp), nrow(coexp))
  bincoexp <- (coexp < 0.05)
  bincoexp[is.na(bincoexp)] <- 0
  rownames(bincoexp) <- genes
  F1_F2_F3_pairwisecoexpr_matrix <- data.frame(matrix(0, nrow = length(all_genes), ncol = ncol(bincoexp)), row.names = all_genes)
  for (g in genes){
    F1_F2_F3_pairwisecoexpr_matrix[g, ] <- bincoexp[g, ]
  }

  ##################################################
  ## CNA
  ##################################################
  cna <- raw.get_tcga_cna(cancer, all_genes, thresholded = FALSE)
  genes <- rownames(cna)
  cna[is.na(cna)] <- 0
  rownames(cna) <- genes
  data_linear_CNA_alltraintest <- data.frame(matrix(0, nrow = length(all_genes), ncol = ncol(cna)), row.names = all_genes)
  for (g in genes){
    data_linear_CNA_alltraintest[g, ] <- cna[g, ]
  }

  ##################################################
  ## Dependency Scores
  ##################################################
  dependency <- raw.get_ccle_rnai_dependency_scores(cancer, all_genes)
  genes <- rownames(dependency)
  mutations <- raw.get_ccle_rnai_mutations(cancer, genes)
  dependency[is.na(dependency)] <- 0
  # Checking for co-dependency with mutated genes
  rownames(dependency) <- genes
  F1_F2_F3_essentiality_matrix <- matrix(1, length(all_genes), length(all_genes))
  rownames(F1_F2_F3_essentiality_matrix) <- colnames(F1_F2_F3_essentiality_matrix) <- all_genes
  for (i in 1:(length(genes)-1)) {
    r <- foreach(j = (i+1):length(genes)) %do% {
      tmp <- c(1, 1)
      for (k in 1:2) {
        if (k == 1) {
          gene1 <- genes[i]
          gene2 <- genes[j]
        } else {
          gene1 <- genes[j]
          gene2 <- genes[i]
        }
        gene1_mutations <- mutations$Tumor_Sample_Barcode[mutations$Hugo_Symbol == gene1]
        if (length(gene1_mutations) == 0) {
          tmp[k] <- 1
        } else {
          x <- as.numeric(dependency[gene2, gene1_mutations])
          y <- as.numeric(dependency[gene2, !(colnames(dependency) %in% gene1_mutations)])
          tmp[k] <- tryCatch({wilcox.test(x, y)$p.value}, error = function(e) 1)
        }
      }
      return(min(tmp[k]))
    }
    r <- unlist(r)
    F1_F2_F3_essentiality_matrix[i, (i+1):ncol(F1_F2_F3_essentiality_matrix)] <- F1_F2_F3_essentiality_matrix[(i+1):nrow(F1_F2_F3_essentiality_matrix), i] <- r
    if (i%%10 == 0) print(paste(i, "of", length(genes)))
  }

  F1_F2_F3_essentiality_matrix <- (F1_F2_F3_essentiality_matrix < 0.05) * 1

  print(paste("finished matrices creation"))
  
  matrices <- list(data_RNA_expression_alltraintest = data_RNA_expression_alltraintest,
       F1_F2_F3_pairwisecoexpr_matrix = F1_F2_F3_pairwisecoexpr_matrix,
       data_linear_CNA_alltraintest = data_linear_CNA_alltraintest,
       F1_F2_F3_essentiality_matrix = F1_F2_F3_essentiality_matrix,
       genes = genes)
  for (m in names(matrices)) {
    write.table(matrices[[m]], paste0(output, m), sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
  }
}

pca_gCMF.predict <- function (train, test, cancer) {
  all_genes <- read_csv(paste0("~/repos/msc-thesis-project/r/data/pca-gCMF/", cancer, "/genes"), 
                        col_names = FALSE)$X1
  output <- paste0("~/repos/msc-thesis-project/r/data/pca-gCMF/", cancer, "/")
  train_and_test <- dplyr::bind_rows(list(train, test))
  
  ##################################################
  ## SL Binary Matrices
  ##################################################
  all_sl <- matrix(NA, ncol = length(all_genes), nrow = length(all_genes))
  rownames(all_sl) <- all_genes
  colnames(all_sl) <- all_genes
  for (i in 1:nrow(train_and_test)) {
    gene1 <- train_and_test$gene1[i]
    gene2 <- train_and_test$gene2[i]
    sl <- ifelse(train_and_test$SL[i] == "Y", 1, 0)
    all_sl[gene1, gene2] <- sl
    all_sl[gene2, gene1] <- sl  
  }
  write.table(all_sl, file = paste0(output, "F1_F2_F3_SL_binary_all"), sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
  
  test_sl <- matrix(NA, ncol = length(all_genes), nrow = length(all_genes))
  rownames(test_sl) <- all_genes
  colnames(test_sl) <- all_genes
  for (i in 1:nrow(test)) {
    gene1 <- test$gene1[i]
    gene2 <- test$gene2[i]
    sl <- ifelse(test$SL[i] == "Y", 1, 0)
    test_sl[gene1, gene2] <- sl
    test_sl[gene2, gene1] <- sl
  }
  write.table(test_sl, paste0(output, "F3_SL_binary_test"), sep = "\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
  
  filename <- paste0(output, "matrix_list_pca.txt")
  fileConn<-file(filename)
  data_files <- c(
    "F3_SL_binary_test	0	1	1",	
    "F1_F2_F3_SL_binary_all	0	1	1",
    "F1_F2_F3_essentiality_matrix	1	1	2",	
    "F1_F2_F3_pairwisecoexpr_matrix	1	1	3",	
    "data_RNA_expression_alltraintest	0	1	4",	
    "data_linear_CNA_alltraintest	1	1	5"	
  )
  writeLines(data_files, fileConn)
  close(fileConn)
  
  gcmf_pca_based_prediction_results <- pca_gCMF(filename, output)
  colnames(gcmf_pca_based_prediction_results) <- c("gene1", "gene2", "score")
  gcmf_pca_based_prediction_results[, 1] <- sapply(as.numeric(gcmf_pca_based_prediction_results[, 1]), function(x) all_genes[x])
  gcmf_pca_based_prediction_results[, 2] <- sapply(as.numeric(gcmf_pca_based_prediction_results[, 2]), function(x) all_genes[x])
  data.frame(gcmf_pca_based_prediction_results)
}

pca_gCMF.filter <- function(pca_gCMF.preds, test) {
  apply(test, 1, function(x, d) {
    gene1 <- x[[1]]
    gene2 <- x[[2]]
    cancer <- x[[3]]
    score <- d$score[d$gene1 == gene1 & d$gene2 == gene2]
    score <- as.numeric(as.character(score))
    return(score)
  }, d = pca_gCMF.preds)
}
