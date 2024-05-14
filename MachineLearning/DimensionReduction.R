### ------ Week 2: Dimension Reduction ------
## ------------- Projection --------------
# MIT Course: High Dimnesional Data Analysis
# Sorce [url]: http://genomicsclass.github.io/book/pages/projections.html
# Student: González Bautista Alejandro

# The projection of some subspcae is the place that
# is closest to the original point
# Simple example with N = 2
if (!require(rafalib)){
  install.packages('rafalib')
  library(rafalib)
}
mypar(1, 1) # It's require to install `rafalib`
plot(c(0,4), c(0,4), xlab = 'Dimension 1',
     ylab = 'Dimension 2', type = 'n')
arrows(0,0,2,3,lwd = 3)
text(2, 3, 'Y', pos = 4,cex = 3)

# Proyection in subspace v
mypar(1, 1)
plot(c(0, 4), c(0, 4), xlab = 'Dimension 1',
     ylab = 'Dimension 2', type = 'n')
arrows(0,0,2,3,lwd = 3)
abline(0,0.5,col = 'red',lwd = 3) # if x = 2c and y= c
# then the slope is 0.5 (y = 0.5x)
text(2, 3, 'Y', pos = 4, cex = 3)
y = c(2, 3)
x = c(2, 1)
cc = crossprod(x, y)/crossprod(x)
segments(x[1]*cc, x[2]*cc, y[1], y[2], lty = 2)
text(x[1]*cc, x[2]*cc, expression(hat(Y)), pos = 4, cex = 3)

# ---------- Rotation ----------------
library(rafalib)
mypar()
plot(c(-2, 4), c(-2, 4), xlab = 'Dimension 1',
     ylab = 'Dimension 2', type = 'n',
     xaxt = 'n', yaxt = 'n', bty = 'n')
text(rep(0, 6), c(c(-2, -1), c(1:4)),
     as.character(c(c(-2, -1),c(1:4))), pos = 2)
text(c(c(-2,-1), c(1:4)), rep(0, 6),
     as.character(c(c(-2, -1),c(1:4))), pos = 1)
abline(v = 0, h = 0)
arrows(0,0,2,3,lwd = 3)
segments(2,0,2,3,lty = 2)
segments(2,3,0,3,lty = 2)
text(2,3,'Y',pos = 4, cex = 3)

# Graphically, we can see that the cordinates
# are the projections to the spaces defined by
# the new basis:
library(rafalib)
mypar()
plot(c(-2,4),c(-2,4),xlab = 'Dimension 1',
     ylab = 'Dimension 2', type = 'n', xaxt = 'n',
     yaxt = 'n', bty = 'n')
text(rep(0,6),c(c(-2,-1),c(1:4)),
     as.character(c(c(-2,-1),c(1:4))), pos = 2)
text(c(c(-2,-1),c(1:4)),rep(0,6),
     as.character(c(c(-2,-1),c(1:4))), pos = 1)
abline(v = 0, h = 0)
abline(0, 1, col = 'red')
abline(0, -1, col = 'red')
arrows(0,0,2,3,lwd = 3)
y <- c(2, 3)
x1 <- c(1, 1) # New basis (i)
x2 <- c(0.5, -0.5) # New basis (j)
c1 <- crossprod(x1, y)/crossprod(x1)
c2 <- crossprod(x2, y)/crossprod(x2)
segments(x1[1]*c1,x1[2]*c1,y[1],y[2],lty = 2)
segments(x2[1]*c2,x2[2]*c2,y[1],y[2],lty = 2)
text(2,3,'Y',pos = 4, cex = 3)

# ---------- Singular Value Decomposition ----------
# Approximate the distance between two dimensional points
# with just one dimension.
# SVD: Y (m x n) = U (m x p) D (p x p) t(V(p x n))
# U -> Orthogonal matrix: Rotation of our data
# V -> Orthogonal matrix
# D -> Diagonal matrix
# p <- min(m, n)

# The variablility (sum of squares to be precise)
# of the columns of t(U)Y = VD are decreasing

library(rafalib)
library(MASS)
n <- 100
y <- t(mvrnorm(n, c(0, 0),
               matrix(c(1,0.95,0.95,1),2,2)))
s <- svd(y)

# Visualize the data of Y
plot(y[1,], y[2,])
# U matrix (2 x 2)
round(sqrt(2) * s$u, 3)

# The plot we showed after the rotation, was showing 
# what we call the principal components
# Columns of the Rotation U^{T}Y
PC1 <- s$d[1] * s$v[,1]
PC2 <- s$d[2] * s$v[,2]
plot(PC1, PC2, xlim = c(-3,3), ylim = c(-3,3))

# How is this useful?
# We will greatky reduce the dimension of V and 
# still be able to reconstruct Y
library(tissuesGeneExpression)
data(tissuesGeneExpression)
set.seed(1)
ind <- sample(nrow(e), 500)
Y <- t(apply(e[ind,], 1, scale)) # Standardize data for illustration

# The `svd` command returns the three matrices
s <- svd(Y)
U <- s$u
V <- s$v
D <- diag(s$d)

# First note that we can in fact recontruct y:
Yhat <- U %*% D %*% t(V)
resid <- Y - Yhat
max(abs(resid))

# If we look at the sum of squares of UD
# We see that the last few are quite close to 0
plot(s$d)

# This implies that the last columns of `V`
# have a very small effect on the reconstruction of `Y`
# SVD is created:
# - - - the columns of V, have less and less influence
# on the reconstruction of Y ("explaining less variance")
k <- ncol(U) - 4
Yhat <- U[,1:k] %*% D[1:k,1:k] %*% t(V[,1:k])
resid <- Yhat - Y
max(abs(resid)) # Yhat ~ Y

# By looking at `d`
# We can obtain a good approximation keeping only 94 columns
# The following plots are useful for seeing
# how much of the variablity is explained by each column
plot(s$d^2/sum(s$d^2)*100,
     ylab = 'Percent variablity explained')

# We can also make a cumulative plot:
plot(cumsum(s$d^2 / sum(s$d^2)) * 100,
     ylab = 'Percent variability explained',
     ylim = c(0, 100), type = 'l')

# Although we start with 189 dimensions, we can approximate
# Y with just 95
k <- 95 # out a possible 189
Yhat <- U[,1:k] %*% D[1:k,1:k] %*% t(V[,1:k])
resid <- Yhat - Y
max(abs(resid))
boxplot(resid, ylim = quantile(Y,c(0.01, 0.99), range = 0))

# Therefore, by using only half as many dimensions
# we retain most of the variablity in our data
1 - var(as.vector(resid))/var(as.vector(Y))
# We say that we explain 96% of the variablity

# Note that we can compute this proportion from D
1 - sum(s$d[1:k]^2)/sum(s$d^2)

## High correlated data
m <- 100
n <- 2
x <- rnorm(m)
e <- rnorm(n*m, 0, 0.01)
Y <- cbind(x,x) + e
cor(Y)

# In this case, the second column adds very little 
# "information" since all entries of Y[,1] - Y[,2]
# are close to 0

# Reporting `rwoMeans(Y)` is even more efficient
# since `Y[,1] - rowMeans(Y)`
# and `Y[,2] - rowMeans(Y)` are even closer to 0.

# The SVD helps us notice that we explain almost all
# the variablity with just this first column
d <- svd(Y)$d
d[1]^2/sum(d^2)

# In cases with many correlated columns, we can achieve great
# dimension reduction
m <- 100
n <- 25
x <- rnorm(m)
e <- rnorm(n*m, 0, 0.01)
Y <- replicate(n, x) + e
d <- svd(Y)$d
cor(Y)
d[1]^2/sum(d^2)


# --------------- Projection Exercices ------------
# We will continue to use this dataset:
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("Biobase")
library(Biobase)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

# Projections Exercises #1
# Suppose you want to make an MA plot
# of the first two samples
# `y = geneExpression`. Which of the following 
# projections of y gives us new coordinates such
# that column 2 vs column 1 is an MA plot?
y <- geneExpression[,1:2]
dim(y)
plot(y[,1],y[,2])
y1 <- y %*% (1/sqrt(2) * matrix(c(1,1,1,-1),2,2))
plot(y1[,1],y1[,2]) # MA Plot
y2 <- y %*% matrix(c(1,1,1,-1),2,2)
plot(y2[,1], y2[,2]) # MA Plot
y3 <- matrix(c(1,1,1,-1),2,2) %*% t(y)
plot(y3[,1], y3[,2])

par(mfrow = c(1,2))
plot(y1[,1],y1[,2]) # MA Plot
plot(y2[,1], y2[,2]) # MA Plot

# SVD Exercises
library(tissuesGeneExpression)
data(tissuesGeneExpression)
# Important Note: When using the SVD in practice
# it is important to note that the solution to SVD 
# is not unique
# U %*% D %*% t(V) = (-U) %*% D %*% t(-V)
s = svd(e)
signflips = sample(c(-1,1), ncol(e), replace = TRUE)
signflips

# Now we switch the sign of each column and check
# that we get the same answer
# `sweep()`

# If x is a matrix and `a` is a vector, then
# `sweep(x,1,y,FUN = '*')` 
  # applies the function `FUN` to each row i
  # `FUN(x[i,],a[i])`, in this case `x[i,] * a[i]`
s$u[1:10,1:3]
length(signflips)
newu <- sweep(s$u, 2, signflips, FUN = '*')
newu[1:10,1:3]
newv <- sweep(s$v, 2, signflips, FUN = '*')
all.equal(s$u %*% diag(s$d) %*% t(s$v),
                       newu %*% diag(s$d) %*% newv)

 
# This is important to know because different implementations
# of the SVD algorithm may give different
# signs, which can lead to the same code resulting
# in different answers when run in different computer systems.

# Compute the SVD of `e`:
m <- rowMeans(e)
U <- s$u
D <- diag(s$d)
V <- s$v
# What is the correlation between the first column
# of U and m?
cor(U[,1],m)

# We saw how the first column relates to the mean of the rows of `e`
# Note that if we change these means, the distances between
# columns do not change.
newmeans <- rnorm(nrow(e)) # Random values we will add to create new means
newe <- e + newmeans # We change the means
sqrt(crossprod(e[,3] - e[,45]))
sqrt(crossprod(newe[,3] - newe[,45]))

# So we might as well make the mean of each row 0
# since it does no help us approximate the column distances.
y <- e - rowMeans(e)
s <- svd(y)

resid <- y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

x <- matrix(rep(c(1,2), each = 5),5,2)
x*c(1:5) # x[i,] * a[i]
sweep(x,1,1:5,'*')

# Which of the following gives us the same 
# as diag(s$d)%*%t(s$v)?
prod1 <- diag(s$d) %*% t(s$v)
prod2 <- s$d * t(s$v)
all.equal(prod1, prod2)
prod3 <- sweep(t(s$v),1,s$d,'*')
all.equal(prod1,prod3)

# If we define vd = t(s$d * t(s$v)), then 
# which of the following is not the same 
# as U %*% D %*% t(V) :
vd = t(s$d * t(s$v))
vd_equal <- s$v * s$d
all.equal(vd,vd_equal)
# tcrossprod(s$u,vd)
eq1 <- s$u %*% (s$d * t(s$v))
eq2 <- s$u %*% t(vd)
eq3 <- tcrossprod(s$u, vd)

# Let z = s$d * t(s$v). 
# We showed a derivation demonstrating that 
# because  is orthogonal, the distance between 
# e[,3] and e[,45] is the same as the distance
# between y[,3] and y[,45], which is the same 
# as z[,3] and z[,45]:
dim(e)
dim(e - rowMeans(e))
z <- s$d * t(s$v)
dim(z)
sqrt(crossprod(e[,3] - e[,45]))
sqrt(crossprod(y[,3] - y[,45]))
sqrt(crossprod(z[,3] - z[,45]))

actual_dist <- sqrt(crossprod(e[,3] - e[,45]))
approx_dist <- sqrt(crossprod(z[1:2,3] - z[1:2,45]))
abs(actual_dist - approx_dist)

# What is the minimum number of dimensions
# we need to use for the approximation in SVD
# to be within 10% or less
actual_dist <- sqrt(crossprod(e[,3] - e[,45]))
approx_dist <- sqrt(crossprod(z[1:7,3] - z[1:7,45]))
abs(actual_dist - approx_dist)/actual_dist * 100

# Compute distances between sample 3 and 
# all other samples:
distances <- sqrt(apply(e[,-3] - e[,3], 2, crossprod))
distances_pca <- sqrt(apply(z[1:2,-3] - z[1:2,3], 2, crossprod))
cor(distances, distances_pca, method = 'spearman')
