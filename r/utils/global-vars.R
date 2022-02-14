# cancer_types <- c("BRCA", "LUAD")
# labels_source <- "discoversl"

discoversl_sysdata_file <- "/scratch/cseale/diff_exp_analysis/DiscoverSL/R/sysdata.rda"
output_dir <- function() paste0("~/repos/SBSL-modelling-and-analysis/processed_data/", labels_source, "/")

get_labels_source <- function() "isle"