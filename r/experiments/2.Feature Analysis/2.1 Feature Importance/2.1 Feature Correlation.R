project_dir <- "~/repos/msc-thesis-project/"
working_dir <- paste0(project_dir, "r/experiments/2.Feature Analysis/2.1 Feature Importance/artifacts/")
setwd(working_dir)
source(paste0(project_dir, "r/utils/train-model.R"))

cancers <- train.get_cancer_types()
train_dataset <- "combined"

# split dataset
train <- train.balance_cancers(train.get_dataset(train_dataset, cancers))
cormat <- round(cor(train[5:ncol(train)]), 2)

library(ggplot2)
library(reshape2)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Melt the correlation matrix
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) +
  coord_fixed()
ggsave("correlation_heatmap.png", width = 8, height = 8, device = "png")

highly_cor <- which(abs(cormat) > .7, arr.ind = T)
highly_cor[apply(highly_cor, 1, function(x) x[1] != x[2]), ]



