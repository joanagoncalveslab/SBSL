project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.6 Training Curves/artifacts/")
setwd(working_dir)
start_time <- Sys.time()
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/logging.R"))

# load parrallel processing
library(doParallel)
library(foreach)
library(pROC)

get_train_and_test_AUC_L0L2 <- function(X, test) {
  train.X <- data.matrix(X[5:ncol(X)])
  train.y <- X$SL
  test.X <- data.matrix(test[5:ncol(test)])
  test.y <- test$SL
  
  m <- train.l0l2(train.X, train.y)
  pred_X <- predict(m, train.X, type = "prob")[, 2]
  X_auc <- auc(roc(train.y, pred_X, direction = "<"))
  pred_test <- predict(m, test.X, type = "prob")[, 2]
  test_auc <- auc(roc(test.y, pred_test, direction = "<"))
  results <- c(X_auc, test_auc)
  names(results) <- c("train", "test")
  results
}

get_train_and_test_AUC_MUVR <- function(X, test) {
  m <- train.MUVR(X)
  pred_X <-  predict(m$Fit$rfFitMax, X, type="prob")[,2]
  X_auc <- auc(roc(X$SL, pred_X, direction = "<"))
  pred_test <- predict(m$Fit$rfFitMax, test, type = "prob")[, 2]
  test_auc <- auc(roc(test$SL, pred_test, direction = "<"))
  results <- c(X_auc, test_auc)
  names(results) <- c("train", "test")
  results
}


generate_learning_curve <- function(train, test, f, filename = NULL) {  
  proportions <- (1:10)/10
  Training_Size <- numeric(length(proportions))
  Test_ROC <- numeric(length(proportions))
  Training_ROC <- numeric(length(proportions))
  Cross_Validation_ROC <- numeric(length(proportions))
  for (i in 1:length(proportions)) {
    p <- proportions[i]
    n <- floor(nrow(train) * p)
    print(paste0("Training on ", p, " of data (", n, " rows)"))
    X <- dplyr::sample_n(train, n)
    # aucs <- get_train_and_test_AUC_L0L2(X, test)
    aucs <- get_train_and_test_AUC_MUVR(X, test)
    print(aucs)
    Training_Size[i] <- n
    Test_ROC[i] <- aucs[["test"]]
    # Cross_Validation_ROC[i] <- mean(m$resample$ROC)
    Training_ROC[i] <- aucs[["train"]]
  }
  
  data <- data.frame(
    Training_Size <- rep(Training_Size, 2), 
    ROC <- c(Test_ROC, Training_ROC),
    Data <- c(rep("Test", length(proportions)), rep("Train", length(proportions)))
  )
   
  # data <- data.frame(
  #   Training_Size <- rep(Training_Size, 3), 
  #   ROC <- c(Test_ROC, Training_ROC, Cross_Validation_ROC),
  #   Data <- c(rep("Test", length(proportions)), rep("Train", length(proportions)), rep("Cross Validation", length(proportions)))
  # )
  
  ggplot(data, aes(x = Training_Size, y = ROC, color = Data)) +
    geom_smooth() + ylim(0.4,1.1) + ggtitle(paste0("Training Curve for ", filename)) + theme(plot.title = element_text(size=10))
  ggsave(paste0(filename, ".png"), width = 6, height = 4)
}

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

cancers <- train.get_cancer_types()

data <- train.balance_cancers(train.get_dataset("combined", cancers))
train_index <- createDataPartition(data$SL, p = .8, 
                                   list = FALSE, 
                                   times = 1)
train <- data[train_index, ]
test <- data[-train_index, ]
preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
f <- train.get_formula(train)
generate_learning_curve(train, test, f, "combined")

for (cancer in cancers) {
    data <- train.balance_cancers(train.get_dataset("combined", cancer))
    train_index <- createDataPartition(data$SL, p = .8, 
                                       list = FALSE, 
                                       times = 1)
    train <- data[train_index, ]
    test <- data[-train_index, ]
    preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
    train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
    test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
    f <- train.get_formula(train)
    generate_learning_curve(train, test, f, cancer)
}

stopCluster(cl)
registerDoSEQ()
log.experiment_run()

end_time <- Sys.time()
print(end_time - start_time)
