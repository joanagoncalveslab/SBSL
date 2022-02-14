project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/0.Dataset Analysis/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/load-labels.R"))

isle <- labels.load("isle")
discoversl <- labels.load("discoversl")

common_pairs_in_isle <- c()
common_pairs_in_discoversl <- c()
for (i in 1:nrow(discoversl)) {
  gene1 <- discoversl$gene1[i]
  gene2 <- discoversl$gene2[i]
  cancer <- discoversl$cancer_type[i]
  
  common_1 <- which(isle$gene1 == gene1 & isle$gene2 == gene2 & isle$cancer_type == cancer)
  common_2 <- which(isle$gene2 == gene1 & isle$gene1 == gene2 & isle$cancer_type == cancer)
  
  if (length(common_1) > 0 | length(common_2) > 0) {
    common_pairs_in_isle <- c(common_pairs_in_isle, common_1) 
    common_pairs_in_isle <- c(common_pairs_in_isle, common_2) 
    common_pairs_in_discoversl <- c(common_pairs_in_discoversl, i)
  }
  if (i%%100==0) print(paste(i, "of", nrow(discoversl)))
}

c_isle <- isle[common_pairs_in_isle, ]
c_discoversl <- discoversl[common_pairs_in_discoversl, ]
print(paste("There are", nrow(c_isle), "pairs in common"))

conf_mat <- table(isle = c_isle$SL, discoverSL = c_discoversl$SL)
print(conf_mat)
#       discoverSL
# isle  0  1
# 0     36 14
# 1     0 19
jaccard_index <- (conf_mat[1,1] + conf_mat[2,2])/sum(conf_mat)
print(paste("The jaccard index is", jaccard_index))
# "The jaccard index is 0.797101449275362"

dis_isle <- c_isle[which(c_isle$SL != c_discoversl$SL), ]
table(dis_isle$cancer_type)
# The disagreement is on 14 BRCA samples
