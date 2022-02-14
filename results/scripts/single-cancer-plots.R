setwd("~/repos/SBSL-modelling-and-analysis/results/")
library(readr)
library(dplyr)
library(xtable)

tested_cancers <- c("BRCA", "LUAD", "OV", "COAD")
AUROCs <- c()
AUPRCs <- c()
AP_3 <- c()
for (cancer in tested_cancers) {
  single <- read_csv(paste0(cancer, "_single.csv"))
  by_model <- single %>% group_by(method)
  AUROCs[[cancer]] <- by_model %>% summarise(
    "AUROC" = sprintf("$%.2f+/-%.2f$", mean(auroc), sd(auroc)),
  )
  AUPRCs[[cancer]] <- by_model %>% summarise(
    "AUPRC" = sprintf("$%.2f+/-%.2f$", mean(auprc), sd(auprc)),
  )
  AP_3[[cancer]] <- by_model %>% summarise(
    "AP_3" = sprintf("$%.2f+/-%.2f$", mean(ap_3), sd(ap_3))
  )
}

AUROC_table <- data.frame(
  Method = AUROCs[["BRCA"]]$method,
  BRCA = AUROCs[["BRCA"]]$AUROC,
  COAD = AUROCs[["COAD"]]$AUROC,
  LUAD = AUROCs[["LUAD"]]$AUROC,
  OV = AUROCs[["OV"]]$AUROC
)
print(xtable(AUROC_table, type = "latex"), include.rownames=FALSE)

AUPRC_table <- data.frame(
  Method = AUPRCs[[cancer]]$method,
  BRCA = AUPRCs[["BRCA"]]$AUPRC,
  COAD = AUPRCs[["COAD"]]$AUPRC,
  LUAD = AUPRCs[["LUAD"]]$AUPRC,
  OV = AUPRCs[["OV"]]$AUPRC
)
print(xtable(AUPRC_table, type = "latex"), include.rownames=FALSE)

AP_3_table <- data.frame(
  Method = AP_3[[cancer]]$method,
  BRCA = AP_3[["BRCA"]]$AP_3,
  COAD = AP_3[["COAD"]]$AP_3,
  LUAD = AP_3[["LUAD"]]$AP_3,
  OV = AP_3[["OV"]]$AP_3
)
print(xtable(AP_3_table, type = "latex"), include.rownames=FALSE)




library(rjson)
library(ggpubr)
source("./scripts/plot_curves.R")

cancers <- c("BRCA", "COAD", "LUAD", "OV")
plots <- list()
tables <- list()
for (cancer in cancers) {
  predictions <- list()
  p <- fromJSON(paste(readLines(paste0("single_", cancer, "_colm_paper.json")), collapse=""))
  print(cancer)
  
  title <- cancer
  filename <- cancer
  V <- 10
  
  plots[[cancer]] <- average_roc_curve_comparisons(p, filename = filename, title = title)
}


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plots$BRCA$boxp)
gg1 <- ggarrange(plots$BRCA$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()), 
                 plots$COAD$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$LUAD$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$OV$rocp + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$BRCA$prcp + theme(legend.position="none", plot.title = element_blank()), 
                 plots$COAD$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$LUAD$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 plots$OV$prcp + theme(legend.position="none", plot.title = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()), 
                 ncol = 4, nrow = 2, legend = "bottom", common.legend = T, widths = c(1.2, 1, 1, 1))
plot(gg1)
ggsave("~/repos/SBSL-modelling-and-analysis/results/figures/combined_roc_prc.pdf", gg1, width = 30, height = 24, units = "cm")

