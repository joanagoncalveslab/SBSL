project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.2 Per Cancer/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
source(paste0(project_dir, "r/utils/train-model.R"))
library(caret)
library(pROC)
library(foreach)
library(ggpubr)
library(xtable)

cancers <- c("BRCA", "COAD", "LUAD", "OV")

plots <- list()
tables <- list()
for (cancer in cancers) {
  predictions <- list()
  for (w in c("single", "baselines")) {
    p <- readRDS(paste0("../", cancer, "_", w, "_predictions.Rdata"))
    predictions$labels <- p$labels
    predictions$preds <- c(predictions$preds, p$preds)
  }
  print(cancer)
  predictions$preds[["RRF"]] <- predictions$preds$`Random Forest`
  predictions$preds$`Random Forest` <- NULL
  for (p in predictions$preds) { print(length(p))}
  
  title <- cancer
  filename <- cancer
  preds <- predictions$preds
  labels <- predictions$labels
  V <- length(preds$`pca-gCMF`)

  plots[[cancer]] <- log.average_roc_curve_comparisons(preds,
                                    labels,
                                    filename = filename,
                                    title = title)
  tables[[cancer]] <- print(log.print_results(preds,
                                              labels,
                                              filename = filename,
                                              title = title))
  # log.logr.regularisation_curves(models$`Elastic Net`,
  #                                filename = filename,
  #                                title = title)
  # print(cancer)
  # print(models$L0L2[[1]]$fit$varnames)
  # for (m in models$L0L2) {
  #   print(plot(m$fit, gamma = m$optimalGamma))
  #   a <- print(m)
  #   print(a[a$lambda == m$optimalLambda & a$gamma == m$optimalGamma, ])
  # }
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plots$BRCA$boxp)
gg1 <- ggarrange(plots$BRCA$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()), 
                 plots$LUAD$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$COAD$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$OV$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$BRCA$prcp + theme(legend.position="none", plot.title = element_blank()), 
                 plots$LUAD$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$COAD$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$OV$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 ncol = 4, nrow = 2, legend = "bottom", common.legend = T, widths = c(1.2, 1, 1, 1))
plot(gg1)
ggsave("combined_roc_prc.pdf", gg1, width = 30, height = 24, units = "cm")


create_latex_tables <- function(filename, tables) { 
  for (d in c("ROC", "PRC", "nCGD")) {
    d_mean <- paste0(d, "_mean")
    d_sd <- paste0(d, "_sd")
    formatted_list <- list()
    for (cancer in cancers) {
      auc_mean <- tables[[cancer]][[d_mean]]
      auc_sd <- tables[[cancer]][[d_sd]]
      r <- paste0("$", round(auc_mean, 2), "\\pm",round(auc_sd, 2), "$")
      names(r) <- tables[[cancer]]$Model
      formatted_list[[cancer]] <- r 
    }
    ggtable <- data.frame(bind_cols(formatted_list))
    rownames(ggtable) <- names(r)
    latex_table <- xtable(ggtable)
    a <- print(xtable(ggtable[order(names(r)), ]), display = c("s", "f", "f", "f", "f"), sanitize.text.function=identity)
    sink(paste0(filename, d, ".txt"))
    cat(str_replace_all(a, "0[.]", "."))
    sink()
  }
}

create_latex_tables("single", tables)

tables <- list()
for (cancer in cancers) {
  predictions <- list()
  for (w in c("all", "unbalanced")) {
    p <- readRDS(paste0("../", cancer, "_", w, "_predictions.Rdata"))
    predictions$labels <- p$labels
    predictions$preds <- c(predictions$preds, p$preds)
  }
  
  title <- cancer
  filename <- cancer
  preds <- predictions$preds
  labels <- predictions$labels
  V <- length(preds$`Elastic Net (All)`)
  
  tables[[cancer]] <- print(log.print_results(preds,
                                              labels,
                                              filename = filename,
                                              title = title))
}
create_latex_tables("all", tables)

create_tsv <- function(tables) {
  for (d in c("ROC", "PRC", "nCGD")) {
    d_mean <- paste0(d, "_mean")
    d_sd <- paste0(d, "_sd")
    formatted_list <- list()
    for (cancer in cancers) {
      auc_mean <- tables[[cancer]][[d_mean]]
      auc_sd <- tables[[cancer]][[d_sd]]
      r <- paste0(round(auc_mean, 2), "\u00B1",round(auc_sd, 2))
      names(r) <- tables[[cancer]]$Model
      formatted_list[[cancer]] <- r 
    }
    ggtable <- data.frame(bind_cols(formatted_list))
    rownames(ggtable) <- names(r)
    write.table(ggtable, file=paste0(d, ".tsv"), quote=FALSE, sep='\t', col.names = NA)
  }
}

tables <- list()
for (cancer in cancers) {
  predictions <- list()
  for (w in c("baselines", "single", "all", "unbalanced")) {
    p <- readRDS(paste0("../", cancer, "_", w, "_predictions.Rdata"))
    predictions$labels <- p$labels
    predictions$preds <- c(predictions$preds, p$preds)
  }
  
  title <- cancer
  filename <- cancer
  preds <- predictions$preds
  labels <- predictions$labels
  V <- length(preds$`Elastic Net (All)`)
  
  tables[[cancer]] <- print(log.print_results(preds,
                                              labels,
                                              filename = filename,
                                              title = title))
}
create_tsv(tables)



model_names <- names(predictions$preds)[c(1,2,3,7)]
models <- list()
for (nm in model_names) {
  models[[nm]] <- readRDS(paste0("../", cancer,"_",nm,".Rdata"))
}



# train models with the selected parameters and measure importance, generate ALE plots, etc

# stopCluster(cl)
# registerDoSEQ()
# for (cancer in cancers) {
#   d <- train.balance_cancers(train.get_dataset("combined", cancer))
#   train_index <- createDataPartition(d$SL, p = .8, 
#                                      list = FALSE, 
#                                      times = 1)
#   train <- d[train_index, ]
#   test <- d[-train_index, ]
#   preProcValues <- caret::preProcess(train[5:ncol(train)], method = c("center", "scale", "nzv"))
#   train <- cbind(train[1:4], predict(preProcValues, train[5:ncol(train)]))
#   test <- cbind(test[1:4], predict(preProcValues, test[5:ncol(test)]))
#   f <- train.get_formula(train)
#   models <- list()
#   models[["L0L2"]] <- train.l0l2(train[5:ncol(train)], train$SL)
#   models[["Elastic Net"]] <- train.logr(train, f)
#   models$MUVR <- train.MUVR(train)
#   log.plot_var_importance_and_ale_plots(models, test, cancer)
# }



