library(readr)
library(readxl)

###########################################
## Generate all the data
###########################################

# source("~/repos/SBSL-modelling-and-analysis/r/features/mutex/discover.R")
# rm(list = ls())

# source("~/repos/SBSL-modelling-and-analysis/r/features/mutex/discoverSL.R")
# rm(list = ls())

# source("~/repos/SBSL-modelling-and-analysis/r/features/mutex/mutex.R")
# rm(list = ls())

# source("~/repos/SBSL-modelling-and-analysis/r/features/mutex/alteration.R")
# rm(list = ls())

# source("~/repos/SBSL-modelling-and-analysis/r/features/pathway_analysis.R")
# rm(list = ls())
# 
# source("~/repos/SBSL-modelling-and-analysis/r/features/generate_deg_features.R")
# rm(list = ls())
# 
# source("~/repos/SBSL-modelling-and-analysis/r/features/expcorr_analysis.R")
# rm(list = ls())
# 
# source("~/repos/SBSL-modelling-and-analysis/r/features/survival_analysis.R")
# rm(list = ls())
# 
# source("~/repos/SBSL-modelling-and-analysis/r/features/dependency_scores_analysis.R")
# rm(list = ls())
# 
# source("~/repos/SBSL-modelling-and-analysis/r/features/GTEx_analysis.R")
# rm(list = ls())


###########################################
## Combine the data
###########################################

source("~/repos/SBSL-modelling-and-analysis/r/utils/load-labels.R")

generate_features_for_dataset <- function(dataset) {
  output_dir <- paste0("~/repos/SBSL-modelling-and-analysis/processed_data/", dataset, "/")
  all_labels <- labels.load(dataset)
  cancer_types <- unique(all_labels$cancer_type)
  # Need to limit cancer types to those which we have generated data for
  cancer_types <- intersect(cancer_types, c("BRCA", "LUAD", "KIRC", "OV"))
  
  results <- list()
  for (C in cancer_types){
    labels <- dplyr::filter(all_labels, all_labels$`cancer_type` == C)
    
    ##########################################################
    # load copathway participation
    ##########################################################
    copathway_participation <- read_delim(paste0(output_dir, C, "_pathway_coparticipation.txt"), 
                                          "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
    
    labels$copathway_participation <- copathway_participation$p.value
    
    ##########################################################
    # mutex values
    ##########################################################
    dsl_mutex <- read_delim(paste0(output_dir, C, "_discoverSL_mutex.txt"), 
                            "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
    labels$dsl_mutex_amp <- dsl_mutex$mutex_amp
    labels$dsl_mutex_del <- dsl_mutex$mutex_del
    labels$dsl_mutex_mut <- dsl_mutex$mutex_mut
    labels$dsl_mutex <- dsl_mutex$mutex
    
    ##########################################################
    # load expression correlation
    ##########################################################
    exp_corr <- read_delim(paste0(output_dir, C, "_exp_corr.txt"), 
                           "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
    labels$exp_corr <- exp_corr$tumour_corr
    labels$exp_corr_pvalue <- exp_corr$tumour_corr.pvalue
    labels$exp_corr_normal <- exp_corr$normal_corr
    labels$exp_corr_normal_pvalue <- exp_corr$normal_corr.pvalue
    
    ##########################################################
    # load differential expressed genes
    ##########################################################
    deg <- read_delim(paste0(output_dir, C, "_deg_pvalues.txt"), 
                      "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
    labels$diff_exp_pvalue <- deg$PValue
    labels$diff_exp_logfc <- deg$logFC
    
    ##########################################################
    # load dependency scores
    ##########################################################
    dependency_scores <- read_delim(paste0(output_dir, C, "_dependency_scores.txt"),
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
    labels$avana_codep <- dependency_scores$avana_codep
    labels$avana_codep_pvalue <- dependency_scores$avana_codep_pvalue
    labels$avana_dep <- dependency_scores$avana_dep
    labels$avana_dep_pvalue <- dependency_scores$avana_dep_pvalue
    labels$avana_avg <- dependency_scores$avana_avg
    labels$d2_codep <- dependency_scores$d2_codep
    labels$d2_codep_pvalue <- dependency_scores$d2_codep_pvalue
    labels$d2_dep <- dependency_scores$d2_dep
    labels$d2_dep_pvalue <- dependency_scores$d2_dep_pvalue
    labels$d2_avg <- dependency_scores$d2_avg
    
    ##########################################################
    # discover mutex scores
    ##########################################################
    discover_mutex <- read_delim(paste0(output_dir, C, "_discover_mutex.txt"),
                                 "\t", escape_double = FALSE,
                                 trim_ws = TRUE)
    labels$discover_mutex <- discover_mutex$discover_mutex
    
    ##########################################################
    # mutex scores
    ##########################################################
    mutex <- read_delim(paste0(output_dir, C, "_mutex.txt"),
                        "\t", escape_double = FALSE,
                        trim_ws = TRUE)
    labels$mutex <- mutex$mutex
    
    ##########################################################
    # survival scores
    ##########################################################
    survival <- read_delim(paste0(output_dir, C, "_survival.txt"),
                           "\t", escape_double = FALSE,
                           trim_ws = TRUE)
    labels$mut_logrank.pval <- survival$mut_logrank_pvals
    labels$mrna_logrank.pval <- survival$mrna_logrank_pvals
    
    ##########################################################
    # survival scores
    ##########################################################
    mutex_alt <- read_delim(paste0(output_dir, C, "_mutex_alt.txt"),
                            "\t", escape_double = FALSE,
                            trim_ws = TRUE)
    labels$mutex_alt <- mutex_alt$mutex_alt
    
    ##########################################################
    # gtex 
    ##########################################################
    gtex <- read_delim(paste0(output_dir, C, "_gtex.txt"),
                       "\t", escape_double = FALSE,
                       trim_ws = TRUE)
    labels$gtex_corr <- gtex$gtex_corr
    labels$gtex_corr.pvalue <- gtex$gtex_corr.pvalue
    
    ##########################################################
    # write completed dataset
    ##########################################################
    results[[C]] <- labels
    print(paste("finished ", C, " features"))
  }
  saveRDS(dplyr::bind_rows(results), file = paste0("~/repos/SBSL-modelling-and-analysis/r/data/", dataset, "_data.Rdata"))
}

generate_features_for_dataset("isle")
generate_features_for_dataset("discoversl")

