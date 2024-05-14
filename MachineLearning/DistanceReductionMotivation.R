# ------------- Distance Reduction --------------
# MIT Course: High Dimnesional Data Analysis
# Sorce [url]: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@2527ed19795c415f8bc5957597175864/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@6a06df24791e45398f9dcefffcd5ce55
# Student: González Bautista Alejandro
n <- 100 # Number of samples
# set.seed(1)

# Multivariate normal distribution
# Covariance matrix
cov_matrix <- matrix(c(1, 0.95, 0.95, 1),2,2)
# Mean vector
mu <- c(0, 0)
if (!require(MASS)){
  install.packages('MASS')
}
library(MASS)
y <- t(mvrnorm(n, mu, cov_matrix))

# Axes tranformation (Matrix A: Rotation)
A <- matrix(c(1/2, 1, 1/2, -1), 2, 2)
z <- A %*% y

# View the transformation
par(mfrow = c(1,2))
plot(y[1,], y[2,])
plot(z[1,], z[2,], xlim = c(-2,2), ylim = c(-3, 4))

d <- dist(t(y)) # 200 x 2 (columns)
d2 <- dist(t(z)) # 200 x 2 (columns)

sd(z[1,])
sd((y[1,] + y[2,])/2)
sd(z[2,])
sd((y[1,] - y[2,]))

# This implies that if we change the transformation above to
A <- 1/sqrt(2) * matrix(c(1, 1, 1 , -1), 2, 2)
solve(A) %*% A # Identity matrix

# This tranformation Preserves Distance
z <- A %*% y
d.y <- dist(t(y)); d.z <- dist(t(z))
as.matrix(d.y)[1, 5]
as.matrix(d.z)[1, 5]

approxd = dist(z[1,]) # Use just one dimension
