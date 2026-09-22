# run_dynamic_cut.R
library(dynamicTreeCut)

my_cutree_function <- function(distM, deep_split = 2, min_size = 5) {
  d <- as.dist(distM)
  hc <- hclust(d, method = "average")
  labels <- cutreeDynamic(
    dendro = hc,
    distM = as.matrix(distM),
    method = "hybrid",
    deepSplit = deep_split,
    minClusterSize = min_size,
    pamRespectsDendro = FALSE
  )
  return(labels)
}