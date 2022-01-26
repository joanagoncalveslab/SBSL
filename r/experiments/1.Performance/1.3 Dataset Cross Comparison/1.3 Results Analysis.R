project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.3 Dataset Cross Comparison/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))
library(caret)

plots <- list()
for (dataset in c("discoversl", "isle")) {
  data <- readRDS(paste0("../", dataset, "_run_data.Rdata"))
  data$preds$DAISY <- NULL
  data$preds$RRF <- data$preds$`Random Forest`
  data$preds$`Random Forest` <- NULL
  title <- ifelse(dataset == "isle", "ISLE/DiscoverSL (BRCA)", "DiscoverSL/ISLE (LUAD)")
  V <- 10
  plots[[dataset]] <- log.average_roc_curve_comparisons(data$preds, data$labels,
                                    filename = dataset,
                                    title = title)
}

plots$isle$rocp

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


mylegend <- g_legend(plots$discoversl$rocp)

gg1 <- ggarrange(plots$isle$rocp, 
                 plots$discoversl$rocp, 
                 boxp,
                 ncol = 3, nrow = 1, widths = c(1, 1, 1), labels = c("A", "B", "C"), common.legend = T, legend= "bottom")
plot(gg1)

ggsave("xdataset_brca.png", plots$isle$rocp, width = 11, height = 11, units = "cm")
ggsave("xdataset_luad.png", plots$discoversl$rocp, width = 11, height = 11, units = "cm")
ggsave("gene_dropout.png", boxp, width = 11, height = 11, units = "cm")

ggsave("xdataset_gene_dropout.pdf", gg1, width = 30, height = 11, units = "cm")

create_tsv <- function(tables) {
  for (d in c("ROC", "PRC", "nCGD")) {
    d_mean <- paste0(d, "_mean")
    d_sd <- paste0(d, "_sd")
    formatted_list <- list()
    for (dataset in c("discoversl", "isle")) {
      auc_mean <- tables[[dataset]][[d_mean]]
      auc_sd <- tables[[dataset]][[d_sd]]
      r <- paste0(round(auc_mean, 2), "\u00B1",round(auc_sd, 2))
      names(r) <- tables[[dataset]]$Model
      formatted_list[[dataset]] <- r 
    }
    ggtable <- data.frame(bind_cols(formatted_list))
    rownames(ggtable) <- names(r)
    write.table(ggtable, file=paste0(d, ".tsv"), quote=FALSE, sep='\t', col.names = NA)
  }
}

tables <- list()
for (dataset in c("discoversl", "isle")) {
  p <- readRDS(paste0("../", dataset, "_run_data.Rdata"))
  title <- ifelse(dataset == "isle", "ISLE/DiscoverSL (BRCA)", "DiscoverSL/ISLE (LUAD)")
  V <- 10
  tables[[dataset]] <- print(log.print_results(p$preds,
                                               p$labels,
                                               filename = filename,
                                               title = title))
}
create_tsv(tables)
