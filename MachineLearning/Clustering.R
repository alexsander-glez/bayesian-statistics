# ---------------------- Clustering ------------------------------
# Basic Machine Learning
# MIT Course: High-Dimensional Data Visualization
# Student: González Bautista Alejandro

# The general idea is to predict or discover
# outcomes from measured predictors.

## Clustering
# The tissue gene expression data
library(tissuesGeneExpression)
data(tissuesGeneExpression)
# View(e)
dim(e) # 189 types of genes
# Let's pretend that we don't know these
# are different tissues and are interested in clustering

# dim(e)[1] Number of samples
# dim(e)[2] Number of genes
# We define the distance between genes
d <- dist(t(e)) # 189 kind of genes

# ------------------- Hierarchical Clustering ----------------------
# With the distance between each pair of samples computed,
# we need clustering algorithms to join them into groups.

# Each sample is assigned to its own group and then the
# algorithms continues iteratively, joining the two most similar
# clusters at each step.

# While we have defined distances between samples, we have not
# yet defined distances between groups
# -> hclust -> -> describes the grouping that were created.

library(rafalib)
mypar()
hc <- hclust(d)
hc # Hierarchical Clustering

plot(hc,labels = tissue, cex = 0.5)
myplclust(hc, labels = tissue, lab.col = as.fumeric(tissue), cex = 0.5)
as.fumeric(tissue)

# Visually, it does seem as if the clustering technique has 
# discovered the tissues.
# - Hierarchical clustering: does not define specific clusters
# - From the dendogram:
#     - We can decipher the distance between any two groups
#     - Looking at the height (at which the two groups split into two).
# To define the clusters:
# - we need to "cut the tree" at some distance
# - group all samples that are within that distance into groups

myplclust(hc, labels = tissue, lab.col = as.fumeric(tissue), cex = 0.5)
abline(h = 120)

hclusters <- cutree(hc, h = 120)
table(true = tissue, cluster = hclusters)

# Find the height that results in the requested number of clusters
hclusters <- cutree(hc, k = 8)
table(true = tissue, cluster = hclusters)

# ** Each tissue is uniquely represented by one of the clusters
# ---------------------- K-means ---------------------------
# We can also cluster with the `kmeans` function to perform k-means clustering.
set.seed(1)
km <- kmeans(t(e[1:2,]), centers = 7) # Gene 1 & 2
e[,1:2] #  GSM11805.CEL.gz GSM11814.CEL.gz
names(km)

mypar(1,2)
plot(e[1,], e[2,], col = as.fumeric(tissue), pch = 16)
plot(e[1,], e[2,], col = km$cluster, pch = 16)
table(true = tissue, cluster = km$cluster)

# MDS plot
km <- kmeans(t(e), centers = 7) # centers: 7 of genes
mds <- cmdscale(d)

mypar(1, 2)
plot(mds[,1], mds[,2])
plot(mds[,1], mds[,2], col = km$cluster, pch = 16)
table(true = tissue, cluster = km$cluster)

# ------------------ Heatmaps ------------------------
library(tissuesGeneExpression)
data(tissuesGeneExpression)
image(e[1:100,]) # First 100 genes
dim(e)
library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:40]

heatmap(e[idx,])

# Change the colors
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
# install.packages('BiocManager')
# BiocManager::install("genefilter")
heatmap(e[idx,], col = hmcol)

if (!require(gplots)){
  install.packages('gplots')
  library(gplots)
}
library(gplots)
cols <- palette(brewer.pal(8, 'Dark2'))[as.fumeric(tissue)]
cbind(colnames(e), tissue, cols)
heatmap.2(e[idx,], labCol = tissue,
          trace = 'none',
          ColSideColors = cols,
          col = hmcol)
