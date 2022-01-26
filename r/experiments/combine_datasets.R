source("~/repos/msc-thesis-project/r/utils/train-model.R")
setwd("~//repos/msc-thesis-project/r/data")

remove_duplicates_within <- function(data) {
  dups <- list()
  agreement <- list()
  for (i in 1:nrow(data)) {
    gene1 <- data$gene1[i]
    gene2 <- data$gene2[i]
    cancer <- data$cancer_type[i]
    SL <- data$SL[i]
    d1 <- which(data$gene1 == gene1 & data$gene2 == gene2 & data$cancer_type == cancer) # this will always return 1 match
    d2 <- which(data$gene2 == gene1 & data$gene1 == gene2 & data$cancer_type == cancer)
    d1 <- d1[d1 != i]
    if (length(d1) > 0 | length(d2) > 0) {
      matches <- unique(c(d1, d2))
      dups[[as.character(i)]] <- matches
      agreement[[as.character(i)]] <- all(data[matches,4] == SL)
    }
    if (i%%100==0) print(paste(i, "of", nrow(data)))
  }
  
  final_data <- data
  if (length(agreement) > 0) {
    agrees <- sum(unlist(agreement))
    disagrees <- sum(!unlist(agreement))
    jaccard_distance <- 1 - (agrees/(agrees+disagrees))
    print("Jaccard distance is :")
    print(jaccard_distance)
    
    removes <- c()
    for (f in names(dups)) {
      f <- as.numeric(f)
      if (!f %in% removes) {
        if (!agreement[[as.character(f)]]) {
          removes <- c(removes, f)
        } 
        removes <- c(removes, dups[[as.character(f)]])
      } 
    }
    
    final_data <- data[-removes, ]
    print(table(final_data$SL, final_data$cancer_type))
  }
  final_data
}

# load data and remove duplicates
isle <- train.get_dataset("isle")
isle <- remove_duplicates_within(isle)
discoversl <- train.get_dataset("discoversl")
discoversl <- remove_duplicates_within(discoversl)
p <- train.duplicate_pairs(isle, discoversl)
isle <- isle[-p, ]
combined <- dplyr::bind_rows(isle, discoversl)
table(combined$SL, combined$cancer_type)
saveRDS(combined, "combined.RData")

