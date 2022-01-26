library(plyr)
library(readr)
library(dplyr)
library(readxl)

source("~/repos/msc-thesis-project/r/utils/load-labels.R")
source("~/repos/msc-thesis-project/r/utils/global-vars.R")

all_labels <- labels.load(labels_source)

for (C in cancer_types) {
  print(paste("Starting", C ,"analysis"))
  results_dir <- paste0("~/repos/msc-thesis-project/processed_data/", C, "_results/")
  diffexp_files <- list.files(results_dir)
  diffexp_files <- diffexp_files[unlist(lapply(diffexp_files, function (x) as.numeric(strsplit(x, "__")[[1]][2]) > 3))]
  
  sample_d <- read_delim(paste0(results_dir, diffexp_files[1]), "\t", escape_double = FALSE, trim_ws = TRUE)
  genes <- sample_d$X1
  cols <- colnames(sample_d)[2:4]
  labels <- dplyr::filter(all_labels, all_labels$cancer_type == C)
  
  diffExp_values <- matrix(NA, dim(labels)[1], 3)
  colnames(diffExp_values) <- cols
  for (i in  1:dim(labels)[1]) {
    diffExp_values[i, ] <- c(0,0,1)
    tryCatch({
      gene1 <- labels[[i, "gene1"]]
      gene2 <- labels[[i, "gene2"]]
      gene1datfile <- diffexp_files[grepl(gene1, diffexp_files)]
      gene1datfile <- read_delim(paste0(results_dir, gene1datfile), "\t", escape_double = FALSE, trim_ws = TRUE)
      result <- gene1datfile[gene1datfile$X1 == gene2, cols]
      if (nrow(result) == 1) {
        diffExp_values[i, ] <- unlist(result)
      } 
    }, error = function(e) {
      diffExp_values[i, ] <- c(0,0,1)
    }, error = function(e) {
      diffExp_values[i, ] <- c(0,0,1)
    })
    
    if (i%%250 == 0) print(paste("Completed", i, "of", dim(labels)[1]))
  }
  write.table(diffExp_values, paste0(output_dir(), C, "_deg_pvalues.txt"), sep = "\t", col.names=cols, row.names = FALSE)
}