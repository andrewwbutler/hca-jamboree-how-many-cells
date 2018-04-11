library(Seurat)
library(rhdf5)
library(RANN)
library(tidyr)
library(dplyr)


readLoom <- function(fpath){
    data <- h5read(fpath,name = "/")
    data_matrix <- t(data$matrix)
    dimnames(data_matrix) <- list(data$row_attrs$gene_names, data$col_attrs$cell_names)
    cluster_ids = data$col_attrs$cluster
    names(cluster_ids) = data$col_attrs$cell_names
    object <- CreateSeuratObject(raw.data = data_matrix)
    object <- NormalizeData(object,display.progress = F,scale.factor = 1e4,normalization.method = "LogNormalize")
    object <- SetIdent(object, names(cluster_ids),cluster_ids)
    return(object)
}


Cosinedist <- function(dataset, cluster1, cluster2){
  # calculate cosine dist
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- AverageExpression(subset, return.seurat=T, show.progress=F)
  res <- proxy::simil(t(as.numeric(subset@data[, 1])),
                      t(as.numeric(subset@data[, 2])),
                      method = "cosine")
  res <- as.vector(res)
  return(res)
}

Spearmandist <- function(dataset, cluster1, cluster2){
  # calculate spearman dist
  # do feature selection first!
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- AverageExpression(subset, return.seurat=T, show.progress=F)

  res <- cor(as.numeric(subset@data[, 1]),
             as.numeric(subset@data[, 2]),
             method = "spearman")
  return(res)
}

Pearsondist <- function(dataset, cluster1, cluster2){
  # calculate pearson dist
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- AverageExpression(subset, return.seurat=T, show.progress=F)

  res <- cor(as.numeric(subset@data[, 1]),
             as.numeric(subset@data[, 2]),
             method = "pearson")
  return(res)
}

CLbalance <- function(dataset, cluster1, cluster2){
  # balance between clusters
  freqtable = table(dataset@ident)
  balance = min(freqtable[cluster1], freqtable[cluster2])/sum(freqtable[c(cluster1, cluster2)])
  return(balance)
}

CLfraction <- function(dataset, fracC1, fracC2){
  # fraction of clusters in dataset
  balance = fracC1 + fracC2
  return(balance)
}

NumCells <- function(dataset, cluster1, cluster2){
  freqtable = table(dataset@ident)
  numCells = freqtable[cluster1] + freqtable[cluster2]
  return(numCells)
}

NumDEgenes <- function(dataset, cluster1, cluster2){
  # library(MAST)
  DEgenes <- abs(rowMeans(dataset@data[,WhichCells(dataset, cluster1)]) - rowMeans(dataset@data[,WhichCells(dataset, cluster2)])) > log(2)
  numDEgenes <- sum(DEgenes)
  return(numDEgenes)
}

NumVARgenes <- function(dataset, cluster1, cluster2){
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- FindVariableGenes(object = subset, mean.function = ExpMean, dispersion.function = LogVMR, do.plot=F)
  numVARgenes <- length(x = subset@var.genes)
  return(numVARgenes)
}

MeanLFC <- function(dataset, cluster1, cluster2){
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- AverageExpression(subset, return.seurat=T, show.progress=F)

  LFC <- abs(dataset@data[,1] - dataset@data[,2])
  LFC <- LFC[!(LFC == Inf | LFC == -Inf | is.na(LFC))]
  return(mean(LFC))
}

VarLFC <- function(dataset, cluster1, cluster2){
  subset <- SubsetData(dataset, cells.use=c(WhichCells(dataset, cluster1), WhichCells(dataset, cluster2)))
  subset <- AverageExpression(subset, return.seurat=T, show.progress=F)

  LFC <- abs(dataset@data[,1] - dataset@data[,2])
  LFC <- LFC[!(LFC == Inf | LFC == -Inf | is.na(LFC))]
  return((sd(LFC))^2)
}


MeanPCDist <- function(object, cluster1, cluster2, pcs=1:30){
    sub.obj <- SubsetData(object, ident.use = c(cluster1, cluster2))
    sub.obj <- FindVariableGenes(sub.obj, do.plot=F, display.progress = F)
    sub.obj@var.genes <- rownames(head(sub.obj@hvg.info, 1000))
    pcs <- 1:min(length(sub.obj@cell.names)-1, length(sub.obj@var.genes)-1, length(pcs))
    # pcs <- 1:min(length(sub.obj@cell.names), length(pcs))
    # sub.obj <- ScaleData(sub.obj, genes.use = sub.obj@var.genes,
    #                   display.progress = FALSE, check.for.norm = FALSE)
    sub.obj <- ScaleData(sub.obj,
                      display.progress = FALSE, check.for.norm = FALSE)
    sub.obj <- RunPCA(sub.obj, pcs.compute = max(pcs), do.print = FALSE)
    c1 <- colMeans(sub.obj@dr$pca@cell.embeddings[WhichCells(sub.obj, cluster1), pcs])
    c2 <- colMeans(sub.obj@dr$pca@cell.embeddings[WhichCells(sub.obj, cluster2), pcs])
    return(as.numeric(dist(rbind(c1, c2))))
}


# ReadDepth <- function(object, cluster1, cluster2){
#     sub.obj <- SubsetData(object, ident.use = c(cluster1, cluster2))
#     return()
# }


# SepDist <- function(object, cluster1 , cluster2, perm.amount = 0.1, n.times = 10) {
#   set.seed = 1
#   orig.dist <- MeanPCDist(object, cluster1, cluster2, pcs = 1:30)
#   new.dist <- c()
#   for (i in 1:n.times) {
#     c1 <- WhichCells(object, cluster1)
#     c2 <- WhichCells(object, cluster2)
#     c1.swap <- sample(c1, length(c1)*perm.amount)
#     c2.swap <- sample(c2, length(c2)*perm.amount)
#     c1 <- setdiff(c1, c1.swap)
#     c1 <- c(c1, c2.swap)
#     c2 <- setdiff(c2,c2.swap)
#     c2 <- c(c2, c1.swap)
#     new.object <- SetIdent(object, cells.use = c1, ident.use = "Cluster1")
#     new.object <- SetIdent(new.object, cells.use = c2, ident.use = "Cluster2")
#     new.dist <- c(new.dist, MeanPCDist(new.object, "Cluster1", "Cluster2", pcs = 1:30))
#   }
#   return(mean(new.dist)-orig.dist)
# }

# MeanPairDistLFC <- function(object, cluster1, cluster2, pcs=1:30){
#   # relative difference of mean pairwise distances within clusters
#   subObj <- SubsetData(object, ident.use = c(cluster1, cluster2))
#   subObj <- FindVariableGenes(subObj, do.plot = F, display.progress = F)
#   subObj@var.genes <- rownames(head(subObj@hvg.info, 1000))
#   subObj <- ScaleData(subObj, genes.use = subObj@var.genes,
#                     display.progress = FALSE, check.for.norm = FALSE)
#   subObj <- RunPCA(subObj, pcs.compute = 50, do.print = FALSE)
#
#   # use
#
#
#   distlist = c()
#   # quick and dirty code: calculate each pair twice
#   cells1 = t(combn(WhichCells(subObj, cluster1), 2))
#   for (i in cells1){
#     for (j in cells1[!(cells1 == i)]){
#       distlist = c(distlist, as.numeric(dist(rbind(subObj@dr$pca@cell.embeddings[i, pcs], subObj@dr$pca@cell.embeddings[j, pcs]))))
#     }
#   }
#
#   meanPairDist1 <- mean(distlist1)
#   meanPairDist2 <- mean(distlist2)
#
#   lfc <- abs(log(meanPairDist1) - log(meanPairDist2))
#
#   return(lfc)
# }


CalcMetrics <- function(dataset, cluster1, fracC1, cluster2, fracC2){
  fracC1 = as.numeric(fracC1)
  fracC2 = as.numeric(fracC2)
  print(paste(cluster1, cluster2, sep=":"))
  metrics = c(NumCells(dataset, cluster1, cluster2),
              Cosinedist(dataset, cluster1, cluster2),
              Spearmandist(dataset, cluster1, cluster2),
              Pearsondist(dataset, cluster1, cluster2),
              CLbalance(dataset, cluster1, cluster2),
              CLfraction(dataset, fracC1, fracC2),
              NumDEgenes(dataset, cluster1, cluster2),
              NumVARgenes(dataset, cluster1, cluster2),
              MeanLFC(dataset, cluster1, cluster2),
              VarLFC(dataset, cluster1, cluster2),
              MeanPCDist(dataset, cluster1, cluster2))
  names(metrics) = c("NumCells", "cosine", "spearman", "pearson", "CLbalance", "CLfraction",
                     "NumDEgenes", "NumVARgenes", "MeanLFC", "VarLFC", "MeanPCDist")
  return(metrics)
}


runLoom <- function(loomFilePath, clusterpairs){
  # takes loomfile and pairs of clusters to calculate their metrics
  dataset <- readLoom(loomFilePath)
  csvFilePath <- gsub("_cells.loom", ".csv", loomFilePath)
  csvFilePath <- gsub("pancreas_", "pancreas_actual_separabilities_", csvFilePath)

  separabilitymatrix <- read.csv(csvFilePath)
  separabilitymatrix <- separabilitymatrix[!is.na(separabilitymatrix$separability),]

  # check whether clusterpairs part of dataset and have separability score and remove missing pair
  delpairs = c()
  for (i in 1:dim(clusterpairs)[1]){
    if (!(clusterpairs$cluster1[i] %in% levels(dataset@ident)) & (clusterpairs$cluster2[i] %in% levels(dataset@ident))){
      delpairs <- c(delpairs, i)
    }
  }

  clp = lapply(separabilitymatrix$clusters, function(x) strsplit(as.vector(x), split=":")[[1]])

  for (i in 1:dim(clusterpairs)[1]){
    if (!(any(sapply(clp, function(x) all(c(as.vector(clusterpairs$cluster1[i]), as.vector(clusterpairs$cluster2[i])) %in% x))))){
      delpairs <- c(delpairs, i)
    }
  }

  delpairs = unique(delpairs)
  print(delpairs)
  if (length(delpairs) > 0){
    clusterpairs = clusterpairs[-delpairs,]
  }

  # check whether file exits or create new one:
  metricfpath = gsub(".loom", "_features.csv", loomFilePath)

  if(!file.exists(metricfpath)){
    metrics = c()
  } else {
    metrics = read.csv(metricfpath)
  }

  separabilitymatrix = read.csv(csvFilePath)
  # shuffle pairs before picking !
  shuffledPairs = sample(1:dim(clusterpairs)[1])
  for (i in shuffledPairs){
    pairmetric <- CalcMetrics(dataset, as.vector(clusterpairs$cluster1[i]),
                                       as.vector(clusterpairs$cluster2[i]))
    pairmetric <- t(t(pairmetric))
    sepb = separabilitymatrix[sapply(clp, function(x) all(c(as.vector(clusterpairs$cluster1[i]), as.vector(clusterpairs$cluster2[i])) %in% x)), "separability"]
    pairmetric <- rbind(pairmetric, sepb)
    rownames(pairmetric)[length(pairmetric)] <- "sepb"
    colnames(pairmetric) <- paste(clusterpairs$cluster1[i], clusterpairs$cluster2[i], sep=":")
    # pairmetric <- sapply(1:dim(clusterpairs)[1], function(i) CalcMetrics(dataset, as.vector(clusterpairs$cluster1[i]),
    #                                                                as.vector(clusterpairs$cluster2[i])))
    metrics <- cbind(metrics, pairmetric)
    write.csv(metrics, file=metricfpath)
  }
  return(metrics)
}



########################
## Calculate Features ##
########################

args = commandArgs(trailingOnly=TRUE)

dataset <- readLoom(args[[1]])

# activated_stellate 0.03 quiescent_stellate 0.02 predictions.tsv

metrics <- CalcMetrics(dataset, args[[2]], args[[3]], args[[4]], args[[5]])

metricsdf = data.frame("NumCells"=seq(1000, 100000, 1000),
                       "cosine"=metrics["cosine"],
                       "spearman"=metrics["spearman"],
                       "pearson"=metrics["pearson"],
                       "CLbalance"=metrics["CLbalance"],
                       "CLfraction"=metrics["CLfraction"],
                       "NumDEgenes"=metrics["NumDEgenes"],
                       "NumVARgenes"=metrics["NumVARgenes"],
                       "MeanLFC"=metrics["MeanLFC"],
                       "VarLFC"=metrics["VarLFC"],
                       "MeanPCDist"=metrics["MeanPCDist"])

metricsdf = t(metricsdf)

write.table(metricsdf, file="features-predicted-separability.csv", sep=",", col.names=F, quote=F)
