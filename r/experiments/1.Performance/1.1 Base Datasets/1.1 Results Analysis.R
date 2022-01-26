project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.1 Base Datasets/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))

library(caret)
library(doParallel)
library(foreach)

data <- readRDS("../run_data.Rdata")
title <- "BRCA, LUAD, OV, and COAD"
V <- 10

log.average_roc_curve_comparisons(data$preds, data$labels,
                                  filename = "All",
                                  title = title)
log.logr.regularisation_curves(data$models$`Elastic Net`,
                               filename = "All",
                               title = title)




# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-25
# As long as, e.g., only continuous predictor variables, as in most gene expression studies, 
# or only variables with the same number of categories are considered in the sample, 
# variable selection with random forest variable importance measures is not affected by our findings. 
# However, in studies where continuous variables, such as the folding energy, 
# are used in combination with categorical information from the neighboring nucleotides, 
# or when categorical predictors, as in amino acid sequence data, vary in their number of 
# categories present in the sample variable selection with random forest variable importance measures 
# is unreliable and may even be misleading. 
# variable importance


# Minimal-optimal (‘Min’) and all-relevant (‘Max’)
# models represent the outer borders of variable selections where the validation performance (in this case,
#                                                                                             number of misclassifications) is minimal. This is in practice determined as having validation performance
# within 5% slack allowance from the actual minimum. The minimal-optimal model thus represents the
# minimal variable set required for optimal method performance, i.e. with the strongest predictors e.g.
# suitable for biomarker discovery. The all-relevant model instead represents the data set with all variables
# with relevant signal-to-noise in relation to the research question: i.e. the strongest predictors and,
# additionally, variables with redundant but not erroneous information. The ‘Mid’ model represents a trade-
#   off between the ‘Min’ and ‘Max’ model and is found at the geometric mean.
# plotVAL(data$models$MUVR[[10]])
# plotStability(data$models$MUVR[[5]], model = "min")
# plotVIP(data$models$MUVR[[5]], model = "max")

# train models with the selected parameters and measure importance, generate ALE plots, etc
stopCluster(cl)
registerDoSEQ()
cancers <- c("COAD")
d <- train.balance_cancers(train.get_dataset("combined", cancers))
train_index <- createDataPartition(d$SL, p = .7, 
                                   list = FALSE, 
                                   times = 1)
train <- d[train_index, ]
test <- d[-train_index, ]
preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
f <- train.get_formula(train)
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

# Comparison of feature importance per model
variables <- as.factor(c(imps$LOL2$feature, imps$`Elastic Net`$feature, imps$MUVR$feature, imps$RRF$feature))
models <- c(rep("L0L2" , nrow(imps$LOL2)), rep("Elastic Net" , nrow(imps$`Elastic Net`)), rep("MUVR" , nrow(imps$MUVR)), rep("RRF" , nrow(imps$RRF)))
value <- c(imps$LOL2$importance, imps$`Elastic Net`$importance, imps$MUVR$importance, imps$RRF$importance)
sd <- c(imps$LOL2$importance.95 - imps$LOL2$importance.05,
        imps$`Elastic Net`$importance.95 - imps$`Elastic Net`$importance.05,
        imps$MUVR$importance.95 - imps$MUVR$importance.05,
        imps$RRF$importance.95 - imps$RRF$importance.05)
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
pL0L2 <- plot_varImp_bar_chart(data1[data1$models == "L0L2", ])
pEN <- plot_varImp_bar_chart(data1[data1$models == "Elastic Net", ])
pMUVR <- plot_varImp_bar_chart(data1[data1$models == "MUVR", ])
pRRF <- plot_varImp_bar_chart(data1[data1$models == "RRF", ])


gg1 <- ggarrange(pL0L2, pEN, pMUVR, pRRF,
                 ncol = 2, nrow = 2)
plot(gg1)

ggsave("COAD_all_var_imp.png", gg1, width = 17, height = 17, units = "cm")

gg2 <- ggarrange(pL0L2, pMUVR,
                 ncol = 2, nrow = 1, widths = c(1, 1))
plot(gg2)
ggsave("COAD_select_var_imp.png", gg2, width = 17, height = 10, units = "cm")

imp_table <- data.frame(
  features = imps$LOL2$feature[order(imps$LOL2$feature)],
  L0L2 = imps$LOL2$importance[order(imps$LOL2$feature)],
  `Elastic Net` = imps$`Elastic Net`$importance[order(imps$`Elastic Net`$feature)],
  `Random Forest` = imps$`Random Forest`$importance[order(imps$`Random Forest`$feature)],
  MUVR = imps$MUVR$importance[order(imps$MUVR$feature)]
)
imp_table[2:5] <- round(imp_table[2:5], 2)
write.table(imp_table, file="imp.tsv", quote=FALSE, sep='\t', col.names = NA)
