library(caret)
library(e1071)
library(readr)
library(ROCR)
library(PRROC)
library(mltools)
library(randomForest)
library(MLmetrics)
library(L0Learn)
library(MUVR)
library(RRF)

"Return the cancer types we will train on"
train.get_cancer_types <- function() c("BRCA", "LUAD", "COAD", "OV")

"Return the metric we use for measuring performance"
train.get_metric <- function() "ROC"

"Return the formula we use for training"
train.get_formula <- function(train, remove_dep_scores = FALSE) {
  vars <- colnames(train)[5:ncol(train)]
  if (remove_dep_scores) {
    vars <- vars[!grepl("(RNAi|CRISPR)", vars)]
  }
  
  f <- as.formula(paste("SL ~", paste(vars, collapse = " + ")))
  f
}

"This function takes a list and returns the combined cancer datasets"
train.get_dataset <- function(dataset=NULL, cancers=NULL) {
  all_data <- readRDS(paste0("~/repos/SBSL-modelling-and-analysis/r/data/", dataset, ".RData"))
  all_data <- na.omit(all_data)
  if (is.null(cancers)) {
    return(all_data)
  }
  cancers <- intersect(cancers, unique(all_data$cancer_type))
  all_data <- all_data[all_data$cancer_type %in% cancers, ]
  if (any(all_data$SL == 1)) {
    all_data$SL <- ifelse(all_data$SL == 1, "Y", "N")
  }
  all_data$SL <- as.factor(all_data$SL)
  
  all_data
}

"This function returns the indices of pairs in d1 which also exist in d2"
train.duplicate_pairs <- function(d1, d2) {
  common_pairs <- c()
  for (i in 1:nrow(d1)) {
    gene1 <- d1$gene1[i]
    gene2 <- d1$gene2[i]
    cancer <- d1$cancer_type[i]
    
    common_1 <- which(d2$gene1 == gene1 & d2$gene2 == gene2 & d2$cancer_type == cancer)
    common_2 <- which(d2$gene2 == gene1 & d2$gene1 == gene2 & d2$cancer_type == cancer)
    
    if (length(common_1) > 0 | length(common_2) > 0) {
      common_pairs <- c(common_pairs, i)
    }
    if (i%%100==0) print(paste(i, "of", nrow(d1)))
  }
  common_pairs
}

train.balance_cancers <- function(data, n = NULL, random_seed = 124) {
  n <- ifelse(is.null(n), min(table(data$cancer_type, data$SL)), n)
  l <- list()
  set.seed(random_seed)
  for (cancer in unique(data$cancer_type)) {
    l[[paste0(cancer, "1")]] <- dplyr::sample_n(data[data$cancer_type == cancer & (data$SL == 1 | data$SL == "Y"), ], n, replace = FALSE)
    l[[paste0(cancer, "0")]] <- dplyr::sample_n(data[data$cancer_type == cancer & (data$SL == 0 | data$SL == "N"), ], n, replace = FALSE)
  }
  dplyr::bind_rows(l)
}

train.balance_classes <- function(data, random_seed = 124) {
  l <- list()
  n <- min(table(data$SL))
  set.seed(random_seed)
  l[["1"]] <- dplyr::sample_n(data[(data$SL == 1 | data$SL == "Y"), ], n, replace = FALSE)
  l[["0"]] <- dplyr::sample_n(data[(data$SL == 0 | data$SL == "N"), ], n, replace = FALSE)
  dplyr::bind_rows(l)
}

"This function creates a logistic regression model using the supplied formula using 10-fold cross validation"
train.logr <- function(train, f){
  train.control <- trainControl(
    method="cv",
    number=10,
    classProbs=TRUE,
    summaryFunction=twoClassSummary,
    verboseIter = TRUE
  )

  train.grid <- expand.grid(alpha = seq(0, 1, length = 20), lambda = seq(0.0001, 1, length = 50))

  m <- caret::train(f
    , data = train
    , method = "glmnet"
    , metric = "ROC"
    , trControl = train.control
    , tuneGrid = train.grid)
  return(m)
}

"This function creates a L0L2 model using the supplied formula using 10-fold cross validation"
train.l0l2 <- function(train, labels) {
  X <- data.matrix(train)
  y <- ifelse(labels == "Y", 1, -1)
  
  cvfit <- L0Learn.cvfit(X, y, 
                         nFolds=10,
                         loss="Logistic",
                         nGamma = 20,
                         nLambda = 50,
                         penalty="L0L2")
  cvMeans <- unlist(lapply(cvfit$cvMeans, min))
  optimalGammaIndex = which.min(cvMeans)
  optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])
  cvfit$optimalLambda = cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
  cvfit$optimalGamma = cvfit$fit$gamma[optimalGammaIndex]
  cvfit
}

"This function creates a random forest model using the supplied formula and 10-fold cross validation"
train.rf <- function(train, f){
  train.control <- trainControl(
    method="cv",
    number=10,
    classProbs=TRUE,
    summaryFunction=twoClassSummary,
    verboseIter = TRUE
  )

  train.grid <- expand.grid(mtry = (4:8))

  m <- caret::train(f
      , data = train
      , method = "rf"
      , tuneLength = 10
      , metric = train.get_metric()
      , trControl = train.control
      , tuneGrid = train.grid)
  return(m)
}

train.RRF <- function(train, f) {
  train.control <- trainControl(
    method="cv",
    number=10,
    classProbs=TRUE,
    summaryFunction=twoClassSummary,
    verboseIter = TRUE
  )

  train.grid <- expand.grid(
    mtry = (4:8),
    coefReg = seq(0.5,1,0.1)
    )
  
  m <- caret::train(f
                , data = train
                , method = "RRFglobal"
                , metric = train.get_metric()
                , trControl = train.control
                , tuneGrid = train.grid)
  return(m)
}

train.MUVR <- function(train) {
  library(MUVR)
  nRep = 5
  nOuter = 10
  varRatio=0.8
  method="RF"
  fitness="AUROC"
  m <- MUVR(X=train[5:ncol(train)],
                     Y=train$SL,
                     nRep = nRep,
                     nOuter=nOuter,
                     fitness=fitness,
                     varRatio = varRatio,
                     method = method)
  m
}

"This function creates a random forest model using the supplied formula and 10-fold cross validation"
train.svm <- function(train, f){
  train.control <- trainControl(
    method="cv",
    number=10,
    classProbs=TRUE,
    summaryFunction=twoClassSummary,
    verboseIter = TRUE
  )
  
  m <- caret::train(f
                    , data = train
                    , method = "svmRadial"
                    , metric = train.get_metric()
                    , trControl = train.control
                    , tuneLength = 10)
  return(m)
}

plot_prc_and_roc <- function(predicitons, test, posClass, negClass) {
  class0.scores <- predicitons
  class0.labels <- test$SL == posClass
  roc <- roc.curve(scores.class0 = class1.scores, weights.class0 = class0.labels, curve = T)
  prc <- pr.curve(scores.class0 = class1.scores, weights.class0 = class0.labels, curve = T)
  return(list(roc, prc))
}

train.analyse <- function(pred, sdev, test, posLabel = "Y", negLabel = "N") {
  curves <- plot_prc_and_roc(pred, test, posLabel, negLabel)
  roc <- curves[[1]]
  prc <- curves[[2]]
  predictions <- data.frame(pred, sdev)
  list(roc = roc, prc = prc, predictions=predictions, test=test$SL)
}


