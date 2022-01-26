setwd("~")
source("~/repos/msc-thesis-project/r/experiments/generate_features.R")
source("~/repos/msc-thesis-project/r/utils/train-model.R")

datasets <- c("isle", "discoversl")
cancers <- c("BRCA", "LUAD", "OV")

num_cores <- parallel::detectCores()
registerDoParallel(cores=num_cores)
for (dataset in datasets) {
  pairs <- train.get_dataset(dataset, cancers)
  pairs <- pairs[pairs$diff_exp_logFC == 1, ]
  pairs <- pairs[-which(colnames(pairs) %in% c("diff_exp_logFC", "diff_exp_pvalue"))]
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
      l[[cancer]] <- generate_diff_expression_scores(gene1, genes, cancer)
      l[[cancer]]$cancer_type <- cancer
    }
    print(paste(g, "of", length(gene1s), "complete"))
    dplyr::bind_rows(l)
  }
  new_pairs <- dplyr::bind_rows(r)

  saveRDS(new_pairs, paste0("~/repos/msc-thesis-project/r/data/", dataset, "_repair_diff_exp.Rdata"))
}
