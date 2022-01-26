project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/1.Performance/1.5 Embeddings/artifacts/images")
setwd(working_dir)
source(paste0(project_dir, "r/utils/logging.R"))

for (cancer in c("BRCA", "LUAD", "OV", "COAD")) {
  data <- readRDS(paste0("../", cancer, "_raw_features_run_data.Rdata"))
  title <- cancer
  filename <- cancer
  V <- 10
  
  range01 <- function(x){(x-(min(x)-1))/((max(x)+1)-(min(x)-1))}
  data$preds$`glm` <- lapply(data$preds$`glm`, range01)
  log.average_roc_curve_comparisons(data$preds, lapply(data$labels, function(x) x == "Y"),
                                    filename = filename,
                                    title = title)
}