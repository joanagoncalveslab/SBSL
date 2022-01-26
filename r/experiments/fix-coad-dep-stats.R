source("~/repos/msc-thesis-project/r/experiments/generate_features.R")
source("~/repos/msc-thesis-project/r/utils/train-model.R")

generate_all_survival_scores_for_gene <- function(gene1, genes, cancer, SL) {
  # details
  df <- data.frame(
    gene1,
    gene2 = genes,
    cancer_type = cancer,
    SL
  )

  df7 <- generate_RNAi_depenency_scores(gene1, genes, cancer)
  df8 <- generate_CRISPR_depenency_scores(gene1, genes, cancer)
  
  combined_df <- Reduce(
    function(x, y, ...) merge(x, y, by=c("gene1", "gene2")), 
    list(df, df7, df8)
  )
  combined_df
}

generate_all_survival_scores_for_dataset <- function(dataset, cancer_list) {
  num_cores <- min(c(parallel::detectCores(), 4))
  registerDoParallel(cores=num_cores)
  pairs <- train.get_dataset(dataset)
  pairs <- pairs[pairs$cancer_type %in% cancer_list, ]
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
      l[[cancer]] <- generate_all_survival_scores_for_gene(gene1, genes, cancer, SL)
    }
    print(paste(g, "of", length(gene1s), "complete"))
    dplyr::bind_rows(l)
  }
  dplyr::bind_rows(r)
}

dataset <- "isle"
cancer_list <- c("COAD")
survival_scores <- generate_all_survival_scores_for_dataset(dataset, cancer_list)
saveRDS(survival_scores, paste0("~/repos/msc-thesis-project/r/data/", dataset, "_COAD_dep_stat_fix.Rdata"))


