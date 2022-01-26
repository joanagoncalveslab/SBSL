source("~/repos/msc-thesis-project/r/utils/load-raw.R")
source("~/repos/msc-thesis-project/r/experiments/baseline/source/discoversl.R")

precalculate_discoverSL_score <- function(cancer, gene1, gene2) {
  myexpression <- raw.get_tcga_expression(cancer)
  myexpression <- cbind(Hugo_Symbol = rownames(myexpression), myexpression)
  mymutdata <- raw.get_tcga_mutations(cancer)[c("Hugo_Symbol", "Tumor_Sample_Barcode",  "Variant_Classification")]
  mycnadata <- raw.get_tcga_cna(cancer)
  mycnadata <- cbind(Hugo_Symbol = rownames(mycnadata), mycnadata)
  print("Pathways, Expression, Mutation, CNA loaded")
  print(paste("Starting sl predictions for", gene1, "in", cancer, "against", gene2))
  #Calculate the parameter SharedPathway
  genepath <- calculatePathway(gene1, gene2)
  #Calculate the parameter DiffExp from TCGA data for a chosen cancer
  dg <- customDataDiffExp(gene1, mymutdata, myexpression)
  if (nrow(dg) == 0) {
    dg <- data.frame(
      Gene1 = gene1,
      Gene2 = gene2,
      logFC = 0,
      logCPM = 0,
      PValue = 1
    )
  }
  #Calculate the parameter Exp.correlation from TCGA data for a chosen cancer
  cg <- customDataExpCorrelation(gene1, gene2, myexpression)
  #Calculate the parameter Mutex from TCGA data for a chosen cancer
  mutdata <- customDataMutex(gene1,gene2, mycnadata, mymutdata)
  #Predict potential synthetic lethal interactors for the chosen primary gene in the chosen cancer
  sl <- predictSL(dg,cg,mutdata,genepath)
  print(paste("Finished sl predictions for", gene1, "in", cancer))
  sl
}

discoverSL.predict <- function(test){
  data_dir <- "~/repos/msc-thesis-project/r/data/"
  discoversl_results <- readRDS(paste0(data_dir, "discoversl_results.RData"))
  no_result_default <- 0.000001 # return value if no precalculated results exists
  apply(test, 1, function(x, d) {
    gene1 <- x[[1]]
    gene2 <- x[[2]]
    cancer <- x[[3]]
    matched_rows <- d[d$gene1 == gene1 & d$gene2 == gene2 & d$cancer_type == cancer, ncol(d) - 1]
    if (length(matched_rows) > 0) {
      return(ifelse(is.na(matched_rows[1]), no_result_default, matched_rows[1]))
    }
    return(no_result_default)
  }, d = discoversl_results)
}
