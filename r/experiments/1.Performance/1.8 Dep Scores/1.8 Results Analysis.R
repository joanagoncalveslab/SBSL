project_dir <- "~/repos/SBSL-modelling-and-analysis/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.8 Dep Scores/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
library(doParallel)
library(foreach)
library(caret)

data <- readRDS("../BRCA_run_data.Rdata")
title <- "Dependency scores model comparison"
V <- 10

data$preds$RRF <- data$preds$`Random Forest`
data$preds$`RRF (No Dep)` <- data$preds$`Random Forest (No Dep)`
data$preds$`Random Forest` <- NULL
data$preds$`Random Forest (No Dep)` <- NULL
plots <- log.average_roc_curve_comparisons(data$preds, data$labels,
                                  filename = "dep_scores",
                                  title = title)

brcabp <- plots$boxp  + scale_x_continuous(limits = c(.5,1), breaks=seq(.5,1,.25)) + coord_flip() + theme_Publication() + scale_color_manual(values=paired.cols) + scale_fill_manual(values=paired.cols) +
  guides(color = guide_legend(ncol=2)) + ggtitle("BRCA") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
brcabp
aggregate(brcabp$data[, 1], list(brcabp$data$Model), mean)
ggsave("BRCA_dep_box_plots.png", brcabp, width = 11, height = 11, units = "cm")



data <- readRDS("../LUAD_run_data.Rdata")
title <- "Dependency scores model comparison"
V <- 10

data$preds$RRF <- data$preds$`Random Forest`
data$preds$`RRF (No Dep)` <- data$preds$`Random Forest (No Dep)`
data$preds$`Random Forest` <- NULL
data$preds$`Random Forest (No Dep)` <- NULL
plots <- log.average_roc_curve_comparisons(data$preds, data$labels,
                                           filename = "dep_scores",
                                           title = title)

luadbp <- plots$boxp+ scale_x_continuous(limits = c(.5,1), breaks=seq(.5,1,.25)) + coord_flip()  + theme_Publication() + scale_color_manual(values=paired.cols) + scale_fill_manual(values=paired.cols) +
  guides(color = guide_legend(ncol=2)) + ggtitle("LUAD") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
luadbp 
aggregate(luadbp$data[, 1], list(luadbp$data$Model), mean)
ggsave("LUAD_dep_box_plots.png", luadbp, width = 11, height = 11, units = "cm")


library(ggpubr)

combp <- ggarrange(brcabp, luadbp, legend = "bottom", common.legend = T, nrow = 2, ncol = 1)
combp 
ggsave("comb_dep_box_plots.pdf", combp, width = 11, height = 11, units = "cm")



create_tsv <- function(tables) {
  for (d in c("ROC", "PRC", "nCGD")) {
    d_mean <- paste0(d, "_mean")
    d_sd <- paste0(d, "_sd")
    formatted_list <- list()
    auc_mean <- tables[[d_mean]]
    auc_sd <- tables[[d_sd]]
    r <- paste0(round(auc_mean, 2), "\u00B1",round(auc_sd, 2))
    names(r) <- tables$Model
    formatted_list <- r 
    write.table(formatted_list, file=paste0("BRCA_", d, ".tsv"), quote=FALSE, sep='\t', col.names = NA)
  }
}


data <- readRDS("../BRCA_run_data.Rdata")
preds <- data$preds
labels <- data$labels
V <- 10
tables <- print(log.print_results(preds,
                                  labels,
                                  filename = filename,
                                  title = title))

create_tsv(tables)





# train models with the selected parameters and measure importance, generate ALE plots, etc
stopCluster(cl)
registerDoSEQ()
cancers <- c("LUAD")
data <- train.balance_cancers(train.get_dataset("combined", cancers))
data$MUTEX <- NULL
vars <- colnames(data)
vars_no_dep <- vars[!grepl("(RNAi|CRISPR)", vars)]
cols_no_dep <- which(colnames(data) %in% vars_no_dep)
train_index <- createDataPartition(data$SL, p = .8,
                                   list = FALSE,
                                   times = 1)
# gene dependency scores
# train <- data[train_index, cols_no_dep]
# test <- data[-train_index, cols_no_dep]
train <- data[train_index, ]
test <- data[-train_index, ]
preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
# remove gene dependency scores
train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
f <- train.get_formula(train)
#glmnet
models <- list()
models[["L0L2"]] <- train.l0l2(train[5:ncol(train)], train$SL)
models[["Elastic Net"]] <- train.logr(train, f)
models$MUVR <- train.MUVR(train)
models[["Random Forest"]] <- train.RRF(train, f)

imps <- list()
imps[["LOL2"]] <- log.plot_var_importance_and_ale_plots(models[["L0L2"]], test, "All", TRUE)
imps[["Elastic Net"]] <- log.plot_var_importance_and_ale_plots(models[["Elastic Net"]], test, "All", FALSE)
imps[["RRF"]] <- log.plot_var_importance_and_ale_plots(models[["Random Forest"]], test, "All", FALSE)
imps[["MUVR"]] <- log.plot_var_importance_and_ale_plots(models$MUVR$Fit$rfFitMax, test, "All", FALSE)

variables <- as.factor(c(imps$LOL2$feature, imps$`Elastic Net`$feature, imps$MUVR$feature, imps$RRF$feature))
models <- c(rep("L0L2" , nrow(imps$LOL2)), rep("Elastic Net" , nrow(imps$`Elastic Net`)), rep("MUVR" , nrow(imps$MUVR)), rep("RRF" , nrow(imps$RRF)))
value <- c(imps$LOL2$importance, imps$`Elastic Net`$importance, imps$MUVR$importance, imps$RRF$importance)
sd <- c(imps$LOL2$importance.95 - imps$LOL2$importance.05,
        imps$`Elastic Net`$importance.95 - imps$`Elastic Net`$importance.05,
        imps$MUVR$importance.95 - imps$MUVR$importance.05,
        imps$RRF$importance.95 - imps$RRF$importance.05)
data1 <- data.frame(models,variables,value, sd)
reorder_within <- function(x, by, within, fun = max, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
plot_varImp_bar_chart <- function (data) {
  # Graph
  ggplot(data, aes(color=variables, y=value, x=reorder_within(variables, value, models))) +
    geom_point() +
    geom_hline(yintercept=1, color="grey55", linetype="dashed") +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) + 
    scale_x_reordered() +
    facet_wrap(~models, nrow = 2, ncol = 2, scales = "free_y") +
    xlab("") +
    scale_color_manual(values=log.many_cols) +
    coord_flip() + theme_Publication() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size=10),
      legend.position = "none") 
}
# pL0L2 <- plot_varImp_bar_chart(data1[data1$models == "L0L2", ])
# pEN <- plot_varImp_bar_chart(data1[data1$models == "Elastic Net", ])
# pMUVR <- plot_varImp_bar_chart(data1[data1$models == "MUVR", ])
# pRRF <- plot_varImp_bar_chart(data1[data1$models == "RRF", ])
# 
# 
# gg1 <- ggarrange(pL0L2, pEN, pMUVR, pRRF,
#                  ncol = 2, nrow = 2)
gg1 <- plot_varImp_bar_chart(data1)
plot(gg1)

ggsave("LUAD_all_var_imp.pdf", gg1, width = 20, height = 20, units = "cm")

gg2 <- ggarrange(pL0L2, pMUVR,
                 ncol = 2, nrow = 1, widths = c(1, 1))
plot(gg2)
ggsave("LUAD_dep_scores_select_var_imp.pdf", gg2, width = 17, height = 10, units = "cm")
stopCluster(cl)
registerDoSEQ()












X <- test[5:ncol(test)]
y = as.numeric(test$SL == "Y")
predictor.rf <- Predictor$new(
  model = muvr_model$Fit$rfFitMax, 
  data = X, 
  y = y, 
  predict.function = pred_fun,
  class = 2,
  type = "prob"
)
imp.rf <- FeatureImp$new(predictor.rf, loss = loss_fun, n.repetitions = 50)$results

# Comparison of feature importance per model
variables <- as.factor(c(imp.rf$feature))
models <- c(rep("Random Forest" , nrow(imp.rf)))
value <- c(imp.rf$importance)
sd <- c(imp.rf$importance.95-imp.rf$importance.05)
data1 <- data.frame(models,variables,value, sd)
plot_varImp_bar_chart <- function (data) {
  # Graph
  ggplot(data, aes(color=variables, y=value, x=reorder(variables, value))) + 
    geom_point() +
    geom_hline(yintercept=1, color="grey55", linetype="dashed") + 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) +
    facet_wrap(~models, nrow = 1, ncol = 2) +
    xlab("") +
    scale_color_manual(values=log.many_cols) + 
    coord_flip() + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size=10),
      legend.position = "none")
}
p2 <- plot_varImp_bar_chart(data1[data1$models == "Random Forest", ])
