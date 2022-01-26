project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/2.Feature Analysis/2.2 Feature Selection/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))
source(paste0(project_dir, "r/utils/analyse-model.R"))
source(paste0(project_dir, "r/utils/baseline.R"))
source(paste0(project_dir, "r/utils/logging.R"))

library(ggplot2)

# get data
cancers <- train.get_cancer_types()
train <- train.balance_cancers(train.get_dataset("combined_train", cancers))
test <- train.balance_cancers(train.get_dataset("combined_test", cancers))


# build models
vars <- colnames(train)[5:ncol(train)]
vars <- vars[-which(vars == "MUTEX")]
vars <- vars[-which(vars == "discoversl_mutex")]
n_vars <- length(vars)
rocs <- numeric(length(vars))
retained_vars <- list()
models <- list()
i <- n_vars
while(length(vars) > 1) {
  print(paste(length(vars), "left"))
  m <- train.logr(train, as.formula(paste("as.factor(SL) ~", paste(vars, collapse = " + "))))
  importance <- varImp(m)$importance
  rocs[i] <- mean(m$results$ROC)
  retained_vars[[i]] <- vars 
  models[[i]] <- m
  vars <- vars[-which(importance == min(importance))]
  i <- i - 1
}
rocs <- rocs[-1]
df <- data.frame(value = rocs, vars = 2:n_vars)

ggplot(data=df, aes(x=vars, y=value)) +
  geom_line()+
  geom_point()
ggsave("logr-var-prec.png", width = 3, height = 3, device = "png")


max_i <- which(rocs==max(rocs))
varImp(models[[max_i]])

log.experiment_run()


