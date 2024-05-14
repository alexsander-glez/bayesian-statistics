## ---------- Hierarchichal Clustering Exercises ----------
# MIT Course: High Dimensional Data Analysis
# url: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@6551b34198154249bb29f977be995fb4/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@7d9c437e50a44996b5a0e9a77bda46b1
# Basic Machine Learning: Clustering 
# Student: González Bautista Alejandro

# RANDOM NUMBER DIFFERENCES
RNGkind()
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# |->| Hierarchical Clustering Exercises #1
# Create a random matrix with no correlation 
set.seed(1)
m <- 1e4
n <- 24
x <- matrix(rnorm(m * n), m, n) # nrows <- m; ncols <- n
colnames(x) <- 1:n

# Run hierarchical clustering on this data
# Create a dendrogram
d <- dist(t(x)) # Cluster the columns
hc <- hclust(d)
class(hc)
plot(hc, cex = 1.0, labels = colnames(x))

# |->| Hierarchical Clustering Exercises #2
# Set the seed at 1 with `set.seed(1)` and replicate
# the creation of this matrix 100 times:
set.seed(1)
m <- 1.0e5
n <- 24
for (i in 1:100){
  x <- matrix(rnorm(m * n), m, n)
}
# Then perform hierarchical clustering
# as in the solution to question 2.4.1
# and find the number of clusters if you use
# `cutree()` at height 143.

# Based on the Monte Carlo simulation, what
# is the standard error of this random variable
set.seed(1)
m <- 1.0e4
n <- 24
nc <- replicate(100,{
  x = matrix(rnorm(m * n), m, n)
  hc = hclust(dist(t(x)))
  length(unique(cutree(hc, h = 143)))
})
plot(table(nc)) ## look at the distribution
popsd(nc)

## --------------------------- K-Means ---------------------------------
# |->| K-means Exercises #1
# Run `kmeans()` with 5 centers for the blood RNA data
library(GSE5859Subset)
data(GSE5859Subset)
k <- 5 # Number of centers
# Set the seed to 10
set.seed(10)
kmeans(sampleInfo, centers = k)
sampleInfo$group
sampleInfo$filename


mds <- cmdscale(dist(t(geneExpression)))
plot(mds[,1], mds[,2])
set.seed(10); k <- 5
km <- kmeans(t(geneExpression), centers <- k)
mypar(1,1)
plot(mds, bg = km$cl, pch = 21)
table(sampleInfo$group, km$cluster)
table(sampleInfo$date, km$cluster)
# Looks better if we re-order:
table(sampleInfo$date, km$cluster)[,c(4,1,5,3,2)]

# ----------------------------- Heatmap ----------------------------------
# |->| Heat Maps Exercises #1
library(GSE5859Subset)
data(GSE5859Subset)
# Pick the 25 genes with the highest across sample
# variance
if(!require(matrixStats)){
  install.packages('matrixStats')
  library(matrixStats)
}
?rowMads # We use mads due to a outlier sample

# heatmap.2 from `gplots` it is a bit
# more customized
library(gplots)
ge <- geneExpression

sampleInfo

# Import libraries
library(rafalib) # as.fumeric(...) one label one color
library(gplots) # heatmap.2(...)
library(RColorBrewer) # brever.pal(...) palette of colors
library(matrixStats) # rowMads(...) sd() for each row
ge <- geneExpression
# Make the colors
ncols <- ncol(ge)
brewer.pal(11, 'RdBu')
cols <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))(25) # Colors for each column
gcol <- brewer.pal(3, 'Dark2')
gcol <- gcol[sampleInfo$g + 1] # Make colors for each group

# Make labels from sampleInfo$date
sampleInfo$date
labelscol <- gsub("2005-", "", sampleInfo$date)

# Pick highly variable genes
sds <- rowMads(ge)
length(sds) == dim(ge)[1]
ind <- order(sds, decreasing = TRUE)[1:25] # 25 higher variable genes
length(ge[ind,])

# Make heatmap
heatmap.2(ge[ind,],
          trace = 'none',
          labRow = geneAnnotation$CHR[ind],
          labCol = labelscol,
          col = cols,
          scale = 'row',
          key = TRUE,
          ColSideColors = gcol)
# Response: A group of chrY genes are higher
# in group 0 and appear to drive the clustering.
# Within those clusters there appears to be clustering by month.

# |->| Heat Maps Exercises #2
# Create a large data set of random data that is completely
# independent of `sampleInfo$group` like this:
library(gplots)
library(matrixStats)
library(genefilter)
library(RColorBrewer)

set.seed(17)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
x <- matrix(rnorm(m * n), m, n)
g <- factor(sampleInfo$g)

# Create two heatmaps with these data.
# Show the group `g` either with labels or colors
# 1. Taking the 50 genes with smallest p-values
# obtained with `rowttests`
# 2. Taking the 50 genes with largest standard deviations
# Make the colors
ncols <- colnames(x)
n <- length(ncols) + 1
cols <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))(25)

# 1. Taking the 50 genes with smallest p-values
# obtained with `rowttests`
pval <- rowttests(x, g) # g: for each factor

# 2. Taking the 50 genes with largest standard deviations
# Make the colors
sds <- rowSds(x)
Indexes <- list(t = order(pval$p.value)[1:50], 
                s = order(-sds)[1:50])

# Make the heatmaps
for (ind in Indexes){
  heatmap.2(x[ind,],
            trace = 'none',
            col = cols,
            scale = 'row',
            labCol = g,
            key = FALSE)
}

# There is no relationship between g and x but with 8,793 tests some 
# will appear significant by chance. Selecting genes with the t-test 
# gives us a deceiving result.