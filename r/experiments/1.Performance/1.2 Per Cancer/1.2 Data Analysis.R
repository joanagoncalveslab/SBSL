project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.3 Dataset Cross Comparison/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/load-raw.R"))

library(ggplot2)
library(ggpubr)
library(reshape2)

# analyse for multicollinearity problem
analyse_multicollinearity <- function() {
  cancers <- train.get_cancer_types()
  data <- train.balance_cancers(train.get_dataset("combined", cancers))
  train_index <- createDataPartition(data$SL, p = .8,
                                     list = FALSE,
                                     times = 1)
  train <- data[train_index, ]
  preProcValues <- preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  f <- train.get_formula(train)
  glm.fit <- glm(f, data = train, family = binomial)
  summary(glm.fit)
  vifs <- car::vif(glm.fit)
  vifs[vifs > 5]
}

gen_heatmap <- function(labels, g) {
  m <- matrix("Unknown", nrow = length(g), ncol = length(g));
  colnames(m) <- g
  row.names(m) <- g
  m[as.matrix(labels[1:2])] <- as.character(labels[[4]])
  m[as.matrix(labels[2:1])] <- as.character(labels[[4]])
  melted <- melt(m)
  p1 <- ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
    scale_fill_manual(values = c("white", "blue", "red")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  p1
}

for (cancer in c("BRCA", "LUAD", "COAD", "OV")) {
  d <- train.get_dataset("combined", cancer)
  g <- union(d$gene1, d$gene2)
  plot(gen_heatmap(d, g))
}


## CCLE mutations
plots <- list()
ggcountsdf <- list()
for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
  d <- train.get_dataset("combined", cancer)
  genes <- unique(union(d$gene1, d$gene2))
  tissues <- c(
    BRCA = "breast",
    LUAD = "lung",
    KIRC = "kidney",
    OV = "ovary",
    COAD = "colorectal",
    LAML = "multiple_myeloma"
  )
  crispr <- readRDS("~/repos/msc-thesis-project/r/data/datasets/avana.RData")
  sample_info <- readRDS("~/repos/msc-thesis-project/r/data/datasets/sample_info.RData")
  samples <- colnames(crispr)
  samples <- intersect(samples, sample_info$DepMap_ID[grepl(tissues[cancer], sample_info$lineage)])  
  mutations <- readRDS("~/repos/msc-thesis-project/r/data/datasets/CCLE_mutations.RData")
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  mutations <- mutations[mutations$DepMap_ID %in% samples, ]
  mutation_counts <- table(mutations$Hugo_Symbol)
  # mutation_counts <- mutation_counts[genes]
  # mutation_counts[is.na(mutation_counts)] <- 0
  # no_mut_genes <- genes[!genes %in% names(mutation_counts)]
  # no_mutation_counts <- numeric(length(no_mut_genes))
  # names(no_mutation_counts) <- no_mut_genes
  mutation_counts <- table(unname(apply(d, 1, function(x) max(c(mutation_counts[x[["gene1"]]], mutation_counts[x[["gene2"]]], 0), na.rm = T))))
  
  counts <- c()
  for (n in names(mutation_counts)) counts <- c(counts, (as.numeric(n) * mutation_counts[n]))
  print(sum(counts)/sum(mutation_counts))
  
  
  # ggcounts <- numeric(11)
  # names(ggcounts) <- c(as.character(0:9), "10+")
  # 
  # for (n in names(ggcounts)[1:10]) {
  #   ggcounts[n] <- mutation_counts[n]
  # }
  # ggcounts[is.na(ggcounts)] <- 0
  # ggcounts["10+"] <- sum(mutation_counts) - sum(ggcounts)
  # 
  # ggcountsdf[[cancer]] <- data.frame(number_of_mutations = names(ggcounts), count = as.numeric(ggcounts))
  # ggcountsdf[[cancer]]$number_of_mutations <- factor(ggcountsdf[[cancer]]$number_of_mutations,
  #                                                    levels = ggcountsdf[[cancer]]$number_of_mutations)
  # plots[[cancer]] <- ggplot(ggcountsdf[[cancer]], aes(x = number_of_mutations, y = count)) + geom_col() + ggtitle(cancer) 
}

gg1 <- ggarrange(plots$BRCA,
                 plots$COAD,
                 plots$LUAD,
                 plots$OV,
                 ncol = 2, nrow = 2)
plot(gg1)


## TCGA mutations
plots <- list()
ggcountsdf <- list()
for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
  d <- train.get_dataset("combined", cancer)
  genes <- unique(union(d$gene1, d$gene2))
  tissues <- c(
    BRCA = "breast",
    LUAD = "lung",
    KIRC = "kidney",
    OV = "ovary",
    COAD = "colorectal",
    LAML = "multiple_myeloma"
  )
  cnv <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-all_thresholded.by_genes.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  maf <- read_delim(paste0("~/repos/msc-thesis-project/raw_data/firehose/20160128-", cancer, "-Mutations-AllSamples.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  
  data.mut <- as.matrix(cnv[, -(1:3)]) != 0
  rownames(data.mut) <- cnv$`Gene Symbol`
  colnames(data.mut) <- substr(colnames(data.mut), 1, 15)
  
  mutated_genes <- maf[maf$Variant_Classification != "Silent", c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  mutated_genes <- mutated_genes[!duplicated(mutated_genes), ]
  mutated_genes$Tumor_Sample_Barcode <- substr(mutated_genes$Tumor_Sample_Barcode, 1, 15)
  mutated_genes <- mutated_genes[mutated_genes$Hugo_Symbol %in% cnv$`Gene Symbol`, ]
  mutated_genes <- mutated_genes[mutated_genes$Tumor_Sample_Barcode %in% colnames(data.mut), ]
  mutated_genes <- mutated_genes[mutated_genes$Hugo_Symbol %in% genes, ]
  
  for (i in 1:nrow(mutated_genes)) {
    gene <- mutated_genes$Hugo_Symbol[i]
    sample <- mutated_genes$Tumor_Sample_Barcode[i]
    data.mut[gene, sample] <- TRUE
    if (i%%1000==0) print(paste(i, "of", nrow(mutated_genes)))
  }
  mutation_counts <- rowSums(data.mut)
  mutation_counts <- table(unname(apply(d, 1, function(x) max(c(mutation_counts[x[["gene1"]]], mutation_counts[x[["gene2"]]], 0), na.rm = T))))

  ggcounts <- numeric(51)
  names(ggcounts) <- c(as.character(0:49), "10+")
  
  for (n in names(ggcounts)[1:50]) {
    ggcounts[n] <- mutation_counts[n]
  }
  ggcounts[is.na(ggcounts)] <- 0
  ggcounts["50+"] <- sum(mutation_counts) - sum(ggcounts)
  
  ggcountsdf[[cancer]] <- data.frame(number_of_mutations = names(ggcounts), count = as.numeric(ggcounts))
  ggcountsdf[[cancer]]$number_of_mutations <- factor(ggcountsdf[[cancer]]$number_of_mutations,
                                                     levels = ggcountsdf[[cancer]]$number_of_mutations)
  plots[[cancer]] <- ggplot(ggcountsdf[[cancer]], aes(x = number_of_mutations, y = count)) + geom_col() + ggtitle(cancer) 
}

gg1 <- ggarrange(plots$BRCA,
                 plots$COAD,
                 plots$LUAD,
                 plots$OV,
                 ncol = 2, nrow = 2)
plot(gg1)

cancer <- "COAD"
data <- train.get_dataset("combined", cancer)
data[1:2, 1:4]
data <- data[, 1:4]
d1 <- data[, c(1,4)]
d2 <- data[, c(2,4)]
colnames(d1)[1] <- "gene"
colnames(d2)[1] <- "gene"
d <- dplyr::bind_rows(list(d1, d2))
head(d)
table(d$gene, d$SL)
r <- table(d$gene, d$SL)
r2 <- rowSums(r)
s <- r2 > 3
names <- names(r2)[s]
ddd <- r[names, 2]/(r[names,1] + r[names,2])


# paste0(names, ":", round(ddd, 2), ", ")


