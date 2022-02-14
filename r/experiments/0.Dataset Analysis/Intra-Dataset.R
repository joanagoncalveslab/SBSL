project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/0.Dataset Analysis/")
source("~/repos/SBSL-modelling-and-analysis/r/utils/train-model.R")

# load data and remove duplicates
cancers <- train.get_cancer_types()
c_type <- "isle"
data <- train.get_dataset(c_type, "BRCA")[, 1:4]

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

agrees <- sum(unlist(agreement))/2
disagrees <- sum(!unlist(agreement))/2
jaccard_distance <- 1 - (agrees/(agrees+disagrees))
print("Jaccard distance is :")
print(jaccard_distance)

keeps <- c()
removes <- c()
for (f in names(dups)) {
  f <- as.numeric(f)
  if (!f %in% removes) {
    keeps <- c(keeps, f)
    removes <- c(removes, dups[[as.character(f)]])
  } 
}

final_data <- data[-removes, ]
print(table(final_data$SL, final_data$cancer_type))

