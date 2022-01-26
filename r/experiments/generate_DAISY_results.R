source("~/repos/msc-thesis-project/r/experiments/baseline/daisy.R")
source("~/repos/msc-thesis-project/r/utils/train-model.R")

dataset <- "discoversl"
cancers <- c("BRCA", "LUAD")
discoversl <- train.get_dataset(dataset, cancers)
genes <- union(discoversl$gene1, discoversl$gene2)
dsl_daisy <- DAISY(discoversl)

dataset <- "isle"
cancers <- c("BRCA", "LUAD", "OV", "COAD", "KIRC")
isle <- train.get_dataset(dataset, cancers)
genes <- union(isle$gene1, isle$gene2)
isle_daisy <- DAISY(isle)

daisy <- dplyr::bind_rows(list(dsl_daisy, isle_daisy))

saveRDS(daisy, paste0("~/repos/msc-thesis-project/r/data/DAISY_results.RData"))
paste("Finished")
