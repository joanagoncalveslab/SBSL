project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.1 Base Datasets/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
library(caret)

# analyse for multicollinearity problem
analyse_multicollinearity <- function() {
  cancers <- train.get_cancer_types()
  data <- train.get_dataset("combined", cancers)
  train_index <- createDataPartition(data$SL, p = .8,
                                     list = FALSE,
                                     times = 1)
  train <- data[train_index, ]
  preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
  train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
  train$SL <- as.numeric(train$SL == "Y")
  f <- train.get_formula(train)
  m.fit <- glm(f, data = train, family = binomial)
  summary(m.fit)
  vifs <- car::vif(m.fit)
  list(VIF = vifs, model = m.fit)
}

v <- analyse_multicollinearity()
v$VIF

summary(v$model)

write.table(v$VIF, file="vif.tsv", quote=FALSE, sep='\t', col.names = NA)