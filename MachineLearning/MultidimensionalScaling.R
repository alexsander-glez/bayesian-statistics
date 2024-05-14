## -------- Multidimensional Scaling --------
# MIT Course: High-Dimensional Data Analysis
# Source [url]: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@6e82f1976ed243b0aa14d68136e1559e/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@a000ee0845cc4e6390f9777c6a799a67
# Student: González Bautista Alejandro

# To simplify the illustration we will only
# consider three tissues
library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)
colind <- tissue%in%c('kidney','colon','liver')
mat <- e[,colind]
group <- factor(tissue[colind])
dim(mat)

# Exploratory step: we wish to know if gene expression
# profiles stored in the columns of `mat` show
# more similarity between tissues than across tissues

# Y = U %*% D %*% t(V) ; D = {d1, d2,...,dn}
# d1 + d2 >> d3 + ... + dn
# Y ~= [U1 U2] %*% matrix(c(d1, 0, 0, d2), 2, 2) %*% t([V1, V2])
# Y ~ [U1 U2] %*% matrix(c(d1, 0, 0, d2), 2, 2) %*% c(vi1 vi2)
# Z = matrix(c(d1, 0, 0, d2), 2, 2) %*% c(vi1 vi2)

s <- svd(mat - rowMeans(mat)) # Reduces the total variation
PC1 <- s$d[1] * s$v[,1]
PC2 <- s$d[2] * s$v[,2]
mypar(1,1)
plot(PC1, PC2, pch = 21, bg = as.numeric(group))
legend('bottomright', levels(group),
       col = seq(along = levels(group)),
       pch = 15, cex = 1.5)

# The accuracy of the approximation depends
# on the proportion of variance explained 
# by the first two principla components.
plot(s$d^2 / sum(s$d^2))
sum(s$d[1:2]^2) / sum(s$d^2) * 100

# Although the first two PCs explain over 50% of the
# variablity, there is plenty of information
# that this plot does not show.
# - Distance between points
PC3 <- s$d[3] * s$v[,3]
PC4 <- s$d[4] * s$v[,4]
mypar(1,1)
plot(PC3,PC4,pch = 21,bg = as.numeric(group))
legend('bottomright',levels(group),
       col = seq(along = levels(group)),
                 pch = 15, cex = 1.5)
sum(s$d[3:4]^2)/sum(s$d^2)

# <<< cmdscale >>>
# Function specifically made for MDS plots
# It takes a distance object as an argument and then uses
# principal component analysis to provide the best 
# approximation to this distance that can be obtained
# with "k" dimensions

d <- dist(t(mat))
mds <- cmdscale(d, k = 2)

mypar()
plot(mds[,1], mds[,2], bg = as.numeric(group), pch = 21,
     xlab = 'First dimension', ylab = 'Second dimension')
legend('bottomleft', levels(group),
       col = seq(along = levels(group)), pch = 15)
pairs(mds)

# These two approaches are equivalent up to an arbitrary
# sign change
mypar(1, 2)
for (i in 1:2){
  plot(mds[,i], s$d[i] * s$v[,i], main = paste('PC',i))
  b = ifelse(cor(mds[,i],s$d[i] * s$v[,i]) > 0, 1, -1) # Sign flip
  abline(0, b) ## b is 1 or -1 depending on the arbitrary sign "flip"
}

# 1. Why the arbitrary sign?
# The SVD is not unique, because:
# (-U) %*% D %*% (-t(V)) = U %*% D %*% t(V)

# 2. Why we substract the mean?
# t({(Y_i - mu) - (Y_j - mu)}){Y_i - mu) - (Y_j - mu)}
# = t({Y_i - Y_j}){Y_i - Y_j}
# Because removing the row averages reduces the
# total variation, it can only make the SVD approximation
# better.

# ----------------- Practice --------------
col.ind <- tissue%in%c('kidney','colon','cerebellum','placenta','liver','hippocampus')
levels(as.factor(tissue))
mat2 <- e[,col.ind]
dim(mat2) # Shape: {22215 rows, 174 columns} rows(samples); columns(features)
group <- factor(tissue[col.ind])


# SVD (Reduce dimensionality)
s2 <- svd(mat2 - rowMeans(mat2))
PC1.2 <- s2$d[1] * s2$v[,1]
PC2.2 <- s2$d[2] * s2$v[,2]
sum(s2$d[1:2]^2)/sum(s2$d^2)
mypar(1,1)
plot(PC1.2,PC2.2,pch = 21, bg = as.numeric(group))
legend('bottomleft', levels(group),
       col = seq(along = levels(group)), pch = 15)

PC3.2 <- s2$d[3] * s2$v[,3]
PC4.2 <- s2$d[4] * s2$v[,4]
sum(s2$d[3:4]^2)/sum(s2$d^2)
mypar(1,1)
plot(PC3.2, PC4.2, pch = 21, bg = as.numeric(group))
legend('topright', levels(group),
       col = seq(along = levels(group)), pch = 15)

# <<< cmdscale >>> command
dim(t(mat2))
d2 <- dist(t(mat2)) 
mds2 <- cmdscale(d2, k = 2)
mypar(1,1)
plot(mds2[,1], mds2[,2],pch = 21, bg = as.numeric(group))
legend('topleft', levels(group), 
       col = seq(along = levels(group)), pch = 15) # group : agrupate the features in one variable

# Create a new variable
z <- diag(d[1:2]) %*% t(s2$v[,1:2])
z <- t(z)
plot(z)
mypar(1,1)
plot(z[,1], z[,2], pch = 21, bg = as.numeric(group))
legend('topleft', levels(group),
       col = seq(along = levels(group)), pch = 15)

mypar(1,3)
plot(PC1.2,PC2.2,pch = 21, bg = as.numeric(group))
legend('bottomleft', levels(group),
       col = seq(along = levels(group)), pch = 15)
plot(mds2[,1], mds2[,2],pch = 21, bg = as.numeric(group))
legend('topleft', levels(group), 
       col = seq(along = levels(group)), pch = 15)
plot(z[,1], z[,2], pch = 21, bg = as.numeric(group))
legend('topleft', levels(group),
       col = seq(along = levels(group)), pch = 15)

# ---------------- MDS Execises ---------------
library(tissuesGeneExpression)
data(tissuesGeneExpression)

# In these exercise we will demonstrate the relationship
# between the SVD and the output of `cmdscale()`
y <- e - rowMeans(e) # Reduces the total variation
s <- svd(y)
z <- s$d * t(s$v)
# We can make an MDS plot:
library(rafalib)
ftissue <- factor(tissue)
mypar(1,1)
plot(z[1,], z[2,], col = as.numeric(ftissue))
legend('topleft', levels(ftissue), col = seq_along(ftissue),
       pch = 1)

# Now run the function `cmdscale()` on the original data:
d <- dist(t(e))
mds <- cmdscale(d)

# 1. What is the correlation between the first row of `z`
# and the first column in `mds`? 
cor(z[1,],mds[,1])

# 2. What is the correlation between the second row of `z`
# and the second column of `mds`?
cor(z[2,], mds[,2])

# Note that the MDS plot is not the same:
library(rafalib)
ftissue = factor(tissue)

mypar(1, 2)
plot(z[1,], z[2,], col = as.numeric(ftissue))
legend('topleft', levels(ftissue), col = seq_along(ftissue), pch = 1)
plot(mds[,1], mds[,2], col = as.numeric(ftissue))

# Load the following dataset:
library(GSE5859Subset)
data(GSE5859Subset)
s <- svd(geneExpression - rowMeans(geneExpression))
z <- s$d * t(s$v)
length(sampleInfo$group)
max(cor(t(z), sampleInfo$group))

# Continue working with the `z` calculated from the
# GSE5859Subset data.
# 5. What is this max correlation?
zt <- t(z) # z(, #dimensions)
max(cor(zt, sampleInfo$group))

# 6. Which dimension of `z` has the seconf highest correlation
# with the outcome `SampleInfo$group`?
sort(cor(zt, sampleInfo$group), decreasing = TRUE)
order(cor(zt, sampleInfo$group), decreasing = TRUE)
cor(zt, sampleInfo$group)[6]

# Note these measurements were made during two months:
sampleInfo$date
# We can extract the month this way
month <- format(sampleInfo$date, "%m")
month <- factor(month)
month

# 7. Which dimension of `z` has the highest correlation with 
# the outcome `month`?
sort(cor(zt, as.numeric(month)), decreasing = TRUE)
order(cor(zt, as.numeric(month)), decreasing = TRUE)
