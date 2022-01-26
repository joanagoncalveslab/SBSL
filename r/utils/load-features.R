library(readr)

for (C in cancer_types) {
  f <- paste0(C, "_features")
  assign(f, read_delim(paste0(output_dir(), f,".txt"), 
                       "\t", escape_double = FALSE, col_types = cols(X1 = col_skip()), 
                       trim_ws = TRUE))
}