source("~/repos/msc-thesis-project/r/experiments/baseline/discoversl.R")
source("~/repos/msc-thesis-project/r/utils/train-model.R")
load("~/repos/DiscoverSL/R/sysdata.rda")

args = commandArgs(trailingOnly=TRUE)
cancer = args[1]
g_index = as.numeric(args[2])
dataset <- args[3]
all_labels <- train.get_dataset(dataset)[,1:3]
labels <- dplyr::filter(all_labels, all_labels$cancer_type == cancer)
gene1s <- unique(labels$gene1)
gene1 <- gene1s[g_index]
pairs <- labels[labels$gene1 == gene1, ]
gene2 <- pairs$gene2

sl <- precalculate_discoverSL_score(cancer, gene1, gene2)
saveRDS(sl, paste0("~/repos/msc-thesis-project/r/data/discoverSL/", dataset, "_", cancer, "_", gene1,".RData"))
paste(paste0("Finished ~/repos/msc-thesis-project/r/data/discoverSL/", dataset, "_", cancer, "_", gene1,".RData"))

# combine files
# setwd("./r/data/discoverSL")
# files <- list.files(".")
# l <- list()
# for (i in 1:length(files)) {
#     l[[i]] <- readRDS(files[i])
#     l[[i]]$cancer_type <- strsplit(files[i], "_")[[1]][2]
# }
# r <- dplyr::bind_rows(l)
# colnames(r)[1:2] <- c("gene1", "gene2")
# saveRDS(r, "../discoversl_results_for_isle.RData")