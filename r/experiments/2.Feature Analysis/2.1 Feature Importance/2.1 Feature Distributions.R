project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/2.Feature Analysis/2.1 Feature Importance/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))

cancers <- c("BRCA", "LUAD", "OV", "COAD")
train_dataset <- "combined"

for (cancer in cancers) {
  # split dataset
  train <- train.get_dataset(train_dataset, cancer)
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  print(preProcValues$method)
  vars <- colnames(train)[5:ncol(train)]
  
  for (v in vars[grepl("CRISPR_dep|RNAi_dep", vars)]) {
    p <- ggplot(data=train, aes(train[[v]])) + 
      geom_histogram() + ggtitle(paste(v, cancer))
    plot(p)
  }
}


# find nearZeroVar features
nxv <- nearZeroVar(train, saveMetrics = TRUE)
nxv[nxv$nzv,]

# find linear combinations
findLinearCombos(train[5:ncol(train)])
