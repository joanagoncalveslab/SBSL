library(caret)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(iml)
library(ROCR)

cols <- c("#ac4a51",
          "#60c352",
          "#795bd0",
          "#589831",
          "#c943a9",
          "#acb742",
          "#cd74de",
          "#3c9754",
          "#d9407f",
          "#6dc38b",
          "#d13a40",
          "#3fc1bf",
          "#d86632",
          "#5b83e2",
          "#d59d35",
          "#8950a2",
          "#617e33",
          "#b296db",
          "#837a33",
          "#5667a8",
          "#d3a26c",
          "#53a2d5",
          "#9a5c2a",
          "#37825e",
          "#ad5b94",
          "#e18079",
          "#964267",
          "#e386b3")

# Feature importance using IML R Package
pred_fun <- function(model, newdata) {
  predict(model, newdata, type = "prob")
}

loss_fun <- function(actual, predicted) {
  1 - Metrics::auc(actual, predicted)
}

calc_imp_over_n_models <- function(train_fun, train, test, f) {
  X <- test[5:ncol(test)]
  y = as.numeric(test$SL == "Y" )
  
  features <- labels(terms(f))[order(labels(terms(f)))]
  imp <- NULL
  for (i in 1:10) {
    m.logr <- train_fun(train, f)
    predictor.logr <- Predictor$new(
      model = m.logr, 
      data = X, 
      y = y, 
      predict.function = pred_fun,
      class = 2,
      type = "prob"
    )
    results <- FeatureImp$new(predictor.logr, loss = loss_fun, n.repetitions = 50)$results
    results <- results[order(results$feature), 2]
    imp <- cbind(imp, results)
  }
  sdev <- apply(imp, 1, sd)
  importance <- apply(imp, 1, mean)
  imp <- data.frame(features, importance, sdev)
  imp <- imp[order(-imp$importance),]
  imp
}

plot_var_importance_and_ale_plots <- function(models, test, cancer) {
  X <- test[5:ncol(test)]
  y = as.numeric(test$SL == "Y" )
  
  m.logr <- models$L0L2
  m.rf <- models$MUVR
  
  predictor.logr <- Predictor$new(
    model = m.logr, 
    data = X, 
    y = y, 
    predict.function = pred_fun,
    class = 2,
    type = "prob"
  )
  imp.logr <- FeatureImp$new(predictor.logr, loss = loss_fun, n.repetitions = 50)$results
  predictor.rf <- Predictor$new(
    model = m.rf, 
    data = X, 
    y = y, 
    predict.function = pred_fun,
    class = 2,
    type = "prob"
  )
  imp.rf <- FeatureImp$new(predictor.rf, loss = loss_fun, n.repetitions = 50)$results
  
  # Comparison of feature importance per model
  variables <- as.factor(c(imp.logr$feature, imp.rf$feature))
  models <- c(rep("Logistic Regression" , nrow(imp.logr)) , rep("Random Forest" , nrow(imp.rf)))
  value <- c(imp.logr$importance, imp.rf$importance)
  sd <- c(imp.logr$importance.95-imp.logr$importance.05, imp.rf$importance.95-imp.rf$importance.05)
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
      scale_color_manual(values=cols) + 
      coord_flip() + 
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none")
  }
  p1 <- plot_varImp_bar_chart(data1[data1$models == "Logistic Regression", ])
  p2 <- plot_varImp_bar_chart(data1[data1$models == "Random Forest", ])
  gg1 <- ggarrange(p1, p2,
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)
  plot(gg1)
  ggsave(paste0(cancer, "_model_variable_importance.png"), width = 8, height = 4, device = "png")
  
  
  # Comparison of new variables
  imp.logr2 <- imp.logr[order(imp.logr$feature), ]
  imp.rf2 <- imp.rf[order(imp.rf$feature), ]
  data2 <- data.frame(variables = as.factor(imp.logr2$feature), logr = imp.logr2$importance, rf = imp.rf2$importance)
  row.names(data2) <- data2$variables
  p3 <- ggplot(data2[grepl("mutex", data2$variables), ], aes(x=logr, y=rf)) + 
    geom_point(color = cols[1], size = 3) + 
    geom_text_repel(aes(label = variables), color = cols[1], size = 3) +
    theme(legend.position = "none") + 
    geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed")
  
  p4 <- ggplot(data2[grepl("(RNAi|CRISPR)", data2$variables), ], aes(x=logr, y=rf)) + 
    geom_point(color = cols[2], size = 3) + 
    geom_text_repel(aes(label = variables), color = cols[2], nudge_y = 0.05, nudge_x = 0.05, size = 3) +
    theme(legend.position = "none") + 
    geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  
  p5 <- ggplot(data2[grepl("(gtex|tumour_corr|normal_corr|diff_exp)", data2$variables), ], aes(x=logr, y=rf)) + 
    geom_point(color = cols[3], size = 3) + 
    geom_text_repel(aes(label = variables), color = cols[3], size = 3) +
    theme(legend.position = "none") + 
    geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  
  p6 <- ggplot(data2[grepl("(logrank)", data2$variables), ], aes(x=logr, y=rf)) + 
    geom_point(color = cols[7], size = 3) + 
    geom_text_repel(aes(label = variables), color = cols[7], size = 3) +
    theme(legend.position = "none") + 
    geom_hline(yintercept=1, color="grey55",  linetype="dashed") + geom_vline(xintercept=1, color="grey55",  linetype="dashed") 
  
  gg2 <- ggarrange(p3, p4, p5, p6,
                   labels = c("C", "D", "E", "F"),
                   ncol = 4, nrow = 1)
  plot(gg2)
  ggsave(paste0(cancer, "_new_var_comparison.png"), width = 12, height = 4, device = "png")
  
  # Effect Size vs P-Value
  pvalue_data <- data2[grepl("(gtex|tumour_corr|normal_corr|RNAi|CRISPR|diff_exp)", data2$variables), ]
  pvalue_data$group <- sapply(pvalue_data$variables, function (x) ifelse(grepl("pvalue", x), "P-Value", "Effect Size"))
  scatterPlot <- ggplot(pvalue_data,aes(x=logr, y=rf, color=group)) + 
    geom_point() + geom_rug() +
    scale_color_manual(values = cols[10:11]) + 
    theme(legend.position=c(0.05,0.99), legend.justification=c(0,1))
  
  plot(scatterPlot)
  ggsave(paste0(cancer, "_effect_vs_pval.png"), width = 4, height = 4, device = "png")
  
  ggarrange(gg1, gg2, ncol = 1, nrow = 2)
  ggsave(paste0(cancer, "_all_var_imp.png"), width = 9, height = 6, device = "png")
  
  ale.rf <- FeatureEffects$new(predictor.rf)
  ale.rf$plot()
  ggsave(paste0(cancer, "_rf_feature_effects.png"), width = 12, height = 9, device = "png")
  ale.logr <- FeatureEffects$new(predictor.logr)
  ale.logr$plot()
  ggsave(paste0(cancer, "_logr_feature_effects.png"), width = 12, height = 9, device = "png")
}
