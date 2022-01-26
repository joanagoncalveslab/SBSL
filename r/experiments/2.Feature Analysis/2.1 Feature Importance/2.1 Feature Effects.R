start_time <- Sys.time()
project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/2.Feature Analysis/2.1 Feature Importance/artifacts/")
models_dir <- paste0(project_dir, "r/experiments/1.Performance/1.2 Per Cancer/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))

library(caret)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(iml)
library(ROCR)

# load parrallel processing
library(foreach)

# Feature importance using IML R Package
pred_fun <- function(model, newdata) {
  predict(model, newdata, type = "prob")
}

calc_ale <- function(m, test, f) {
  X <- test[5:ncol(test)]
  X <- X[-which(colnames(X) %in% train.get_excluded_features())]
  y = as.numeric(test$SL == "Y" )
  
  predictor.logr <- Predictor$new(
    model = m, 
    data = X, 
    y = y, 
    predict.function = pred_fun,
    class = 2,
    type = "prob"
  )
  ale <- FeatureEffects$new(predictor.logr)
  ale
}

# Datasets
cancers <- train.get_cancer_types()
cancers <- "BRCA"
train_dataset <- "combined_train"
test_dataset <- "combined_test"

train <- train.balance_cancers(train.get_dataset(train_dataset, cancers))
test <- train.balance_cancers(train.get_dataset(test_dataset, cancers))
f <- train.get_formula(train)

X <- test[5:ncol(test)]
X <- X[-which(colnames(X) %in% train.get_excluded_features())]
y = as.numeric(test$SL == "Y" )

rf <- train.rf(train, f)
predictor.rf <- Predictor$new(
  model = rf, 
  data = X, 
  y = y, 
  predict.function = pred_fun,
  class = 2,
  type = "prob"
)
logr <- train.logr(train, f)
predictor.logr <- Predictor$new(
  model = logr, 
  data = X, 
  y = y, 
  predict.function = pred_fun,
  class = 2,
  type = "prob"
)

ale.rf <- FeatureEffects$new(predictor.rf)
ale.rf$plot()
ale.logr <- FeatureEffects$new(predictor.logr)
ale.logr$plot()

interact.rf <- Interaction$new(predictor.rf)
plot(interact.rf)
interact.logr <- Interaction$new(predictor.logr)
plot(interact.logr)

interact <- Interaction$new(predictor.logr)
plot(interact)

end_time <- Sys.time()
print(end_time - start_time)
