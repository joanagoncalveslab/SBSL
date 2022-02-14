project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.2 Per Cancer/artifacts/")
setwd(working_dir)
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/train-model.R"))

k <- 10
fold_indices <- list()
for (cancer in c("BRCA", "COAD", "LUAD", "OV")){
  data <- train.get_dataset("combined", cancer)
  
  
  for (i in 1:k) {
    fold_indices[[cancer]][[i]] <- createDataPartition(data$SL, p = .7, 
                        list = FALSE, 
                        times = 1)
  }

  for (i in 1:k) {
    f <- fold_indices[[cancer]][[i]]
    d <- data[f, ]
    if (i == k) print(table(d$SL, d$cancer_type))
  }
}
saveRDS(fold_indices, "train_test_indices.Rdata")