source("metrics.R")

# run for pancreas:
sepscores <- read.csv("/home/jovyan/scratch/Group1/downsampled_data/pancreas/pancreas_actual_separabilities_0.95_9.csv")
sepscores <- sepscores[!is.na(sepscores$separability), ]
clusterpairs <- t(sapply(sepscores$clusters, function(clpair) strsplit(as.vector(clpair), ":")[[1]]))
colnames(clusterpairs) <- c("cluster1", "cluster2")
clusterpairs = as.data.frame(clusterpairs)

csvFilePaths <- list.files(path="/home/jovyan/scratch/Group1/downsampled_data/pancreas", pattern="actual_separabilities_*.*csv$", full.names=T)
loomFilePaths <- gsub(".csv", "_cells.loom", csvFilePaths)
loomFilePaths <- gsub("actual_separabilities_", "", loomFilePaths)

parallel::mclapply(loomFilePaths, function(loomFilePath) runLoom(loomFilePath, clusterpairs), mc.cores=16)


# ls -l *features.csv | wc -l
# ls -l *_features.csv | wc -l
# sequencing depth


# run for bipolar:
sepscores <- read.csv("/home/jovyan/scratch/Group1/downsampled_data/bipolar/bipolar_actual_separabilities_0.7_9.csv")
sepscores <- sepscores[!is.na(sepscores$separability), ]
clusterpairs <- t(sapply(sepscores$clusters, function(clpair) strsplit(as.vector(clpair), ":")[[1]]))
colnames(clusterpairs) <- c("cluster1", "cluster2")
clusterpairs = as.data.frame(clusterpairs)

csvFilePaths <- list.files(path="/home/jovyan/scratch/Group1/downsampled_data/bipolar", pattern="actual_separabilities_*.*csv$", full.names=T)
loomFilePaths <- gsub(".csv", "_cells.loom", csvFilePaths)
loomFilePaths <- gsub("actual_separabilities_", "", loomFilePaths)

parallel::mclapply(loomFilePaths, function(loomFilePath) runLoom(loomFilePath, clusterpairs), mc.cores=16)



# run for pbmc:
sepscores <- read.csv("/home/jovyan/scratch/Group1/downsampled_data/pbmc/pbmc_actual_separabilities_0.7_9.csv")
sepscores <- sepscores[!is.na(sepscores$separability), ]
clusterpairs <- t(sapply(sepscores$clusters, function(clpair) strsplit(as.vector(clpair), ":")[[1]]))
colnames(clusterpairs) <- c("cluster1", "cluster2")
clusterpairs = as.data.frame(clusterpairs)

csvFilePaths <- list.files(path="/home/jovyan/scratch/Group1/downsampled_data/pbmc", pattern="actual_separabilities_*.*csv$", full.names=T)
loomFilePaths <- gsub(".csv", "_cells.loom", csvFilePaths)
loomFilePaths <- gsub("actual_separabilities_", "", loomFilePaths)

parallel::mclapply(loomFilePaths, function(loomFilePath) runLoom(loomFilePath, clusterpairs), mc.cores=16)
