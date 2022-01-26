# Different Ideas for dealing with the class imbalance
# https://www.kdnuggets.com/2017/06/7-techniques-handle-imbalanced-data.html

project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.4 Ensemble/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/experiments/baseline/pca-gCMF.R"))
source(paste0(project_dir, "r/experiments/baseline/daisy.R"))
source(paste0(project_dir, "r/experiments/baseline/discoversl.R"))

library(doParallel)
library(foreach)
library(caret)
library(caretEnsemble)

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

cancers <- train.get_cancer_types()
train_dataset <- "combined_train"
test_dataset <- "combined_test"
do.call(file.remove, list(list.files(working_dir, full.names = TRUE)))

# split dataset
train <- train.balance_cancers(train.get_dataset(train_dataset, cancers))
test <- train.balance_cancers(train.get_dataset(test_dataset, cancers))

# Dataset Makeup
log.dataset_makeup(train, filename = train_dataset)
log.dataset_makeup(test, filename = test_dataset)

# get predictors, remove MUTEX because it is a useless feature
vars <- colnames(train)[5:ncol(train)]
vars <- vars[-which(vars == "MUTEX")]
f <- as.formula(paste("as.factor(SL) ~", paste(vars, collapse = " + ")))


# train and save models
# train and save models
train.control <- trainControl(
  method="cv",
  number=10,
  classProbs=TRUE,
  summaryFunction=prSummary,
  verboseIter = TRUE
)

train.grid <- expand.grid(alpha = 0:1, lambda = seq(0.0001, 1, length = 100))

m <- train(f
    , data = train
    , method = "glmnet"
    , metric = "AUC"
    , trControl = train.control
    , tuneGrid = train.grid)

train.analyse(m, test)

stopCluster(cl)
log.experiment_run()
