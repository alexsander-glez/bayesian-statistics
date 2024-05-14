## -------- Principal Component Analysis ------------
# MIT Course: High-Dimensional Data Analysis
# Source [url]: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@6e82f1976ed243b0aa14d68136e1559e/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@2632eb3ada2347fabf81092ac68ee653
# Student: González Bautista Alejandro
if (!require("rafalib")){
  install.packages("rafalib")
  library(rafalib)
}
library(MASS)

n <- 100 # Number of samples
# set.seed(1)

# Multivariate normal distribution
# Covariance matrix
cov_matrix <- matrix(c(1, 0.95, 0.95, 1),2,2)
# Mean vector
mu <- c(0, 0)
Y <- t(mvrnorm(n, mu, cov_matrix))
dim(Y)

# Maximized t({t(u_1)Y}) %*% {t(u_1)Y}
# This can be viewed as a projection of each sample or
# column of Y into the subspace spanned by u_1.

# So we are looking for a transformation in which the coordinates
# show *high variablility*
# Let's try u = t(1, 0)
thelim <- c(-3, 3)

mypar(1,1)
plot(t(Y), xlim = thelim, ylim = thelim,
     main = paste('Sum of Squares: ', round(crossprod(Y[,1]),1)))
abline(h = 0)
apply(Y, 2, function(y) segments(y[1],0,y[1],y[2], lty = 2))
points(Y[1,], rep(0, ncol(Y)), col = 2, pch = 16, cex = 0.65)


# Let's try u = t(1/sqrt(2), -1/sqrt(2)) 
# t(u) %*% u = 1
u <- matrix(c(1, -1)/sqrt(2), ncol = 1)
w <- t(u) %*% Y # Principal component
mypar(1, 1)
plot(t(Y),
     main = paste('Sum of Squares: ',
                  round(tcrossprod(w), 1)),
     xlim = thelim, ylim = thelim)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(0, -1, col = 2)
Z <- u %*% w # Projection of w over u

# w %*% t(w) == tcrossprod(w) == t(t(w) %*% w)
# t(w) %*% w == crossprod(w)
for(i in seq(along = w)){
  segments(Y[1,i], Y[2,i], Z[1,i], Z[2,i], lty = 2)
}
points(t(Z), col = 2, pch = 16, cex = 0.5)

# Let's try u = t(1/sqrt(2), 1/sqrt(2))
u <- matrix(c(1, 1)/sqrt(2), ncol = 1)
t(u) %*% u # Orthogonal matrix
w <- t(u) %*% Y # Projection over u
dim(w)
mypar(1, 1)
plot(t(Y),
     main = paste('Sum of Squares:', round(tcrossprod(w), 1)),
     xlim = thelim, ylim = thelim)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(0, 1, col = 2)
Z <- u %*% w
for(i in seq(along = w)){
  segments(Y[1,i], Y[2,i], Z[1,i], Z[2,i], lty = 2)
}
points(t(Z), col = 2, pch = 16, cex = 0.5)

# The Principal Components
# v {PC} -> maximizes the sum of squares
# The orthogonal vector that maximizes the sum of
# squares: t({t(u_1) %*% Y}) %*% {t(u_1) %*% Y}
# t(u) %*% Y = diag(d) %*% t(v)
# For first PC: t(u_1) %*% Y = diag(d_1) %*% t(v_1)
# To obtain the second PC:
# r = Y - t(u_1) %*% Y %*% v_1

# <<< prcomp >>>
pc <- prcomp(t(Y))
# Produces the same results as the SVD up to arbitrary
# sign flips
s <- svd(Y - rowMeans(Y))
mypar(1,2)
for(i in 1:nrow(Y)){
  plot(pc$x[,i], s$d[i] * s$v[,i])
  print(cor(pc$x[,i], s$d[i] * s$v[,i]))
}

# The loadings can be found this way:
pc$rotation
s$u

# The equivalent of the variance explained 
# is included in the:
pc$sdev