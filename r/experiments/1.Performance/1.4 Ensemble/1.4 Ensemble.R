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

genes <- unique(c(train$gene1, train$gene2, test$gene1, test$gene2))

# Dataset Makeup
log.dataset_makeup(train, filename = train_dataset)
log.dataset_makeup(test, filename = test_dataset)

# get predictors, remove MUTEX because it is a useless feature
vars <- colnames(train)[5:ncol(train)]
vars <- vars[-which(vars == "MUTEX")]
f <- as.formula(paste("as.factor(SL) ~", paste(vars, collapse = " + ")))

model_weights <- ifelse(train$SL == "Y",(1/table(train$SL)[1]) * 0.5, (1/table(train$SL)[2]) * 0.5)

# train and save models
my_control <- trainControl(
  method="boot",
  number=25,
  savePredictions="final",
  classProbs=TRUE,
  index=createResample(train$SL, 25),
  summaryFunction=prSummary,
  search="random"
)

model_list <- caretList(
  f, data=train,
  trControl=my_control,
  metric="AUC",
  methodList=c("glm", "rf", "gbm")
  )

greedy_ensemble <- caretEnsemble(
  model_list, 
  metric="AUC",
  trControl=trainControl(
    number=2,
    summaryFunction=prSummary,
    classProbs=TRUE
    ))



# get results for each model
library("caTools")
model_preds <- lapply(model_list, predict, newdata=test, type="prob")
model_preds <- lapply(model_preds, function(x) x[,"Y"])
model_preds <- data.frame(model_preds)
ens_preds <- predict(greedy_ensemble, newdata=test, type="prob")
model_preds$ensemble <- ens_preds
caTools::colAUC(model_preds, test$SL)

glm_ensemble <- caretStack(
  model_list,
  method="rpart",
  metric="AUC",
  trControl=trainControl(
    method="boot",
    number=10,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=prSummary
  )
) 

model_preds2 <- model_preds
model_preds2$ensemble <- predict(glm_ensemble, newdata=test, type="prob")
CF <- coef(glm_ensemble$ens_model$finalModel)[-1]
colAUC(model_preds2, test$SL)

stopCluster(cl)
registerDoSEQ()
log.experiment_run()


