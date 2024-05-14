# -------------- Distance in High Dimensions -------------
## "High-Dimensional Data Analysis"
# MIT Course: Edx platform
# Class link: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@2527ed19795c415f8bc5957597175864/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@68517e3bbdd7489a829acdd2e954dd38
# Student: González Bautista Alejandro

if (!require(devtools)){
  install.packages('devtools')
  library(devtools)
}

install_github("genomicsclass/tissuesGeneExpression")
# We need to install Rtools for use this command

# The data represent RNA expression levels for **eight tissues**
# each with several individuals
library(tissuesGeneExpression)
data(tissuesGeneExpression)
# View(e)
dim(e)

table(tissue)
sum(table(tissue)) # tissue[i] tells us what tissue is represented by e[,i]

# We are interested in describing distance between samples in the context of
# this dataset (genes that behave similarity across samples)
# Samples 1, 2: kidneys
# Sample 87: colon
x <- e[,1]
y <- e[,2]
z <- e[,87]

# Distance 1: vectorial difference
(dist_xy_1 <- sqrt(sum((x - y)^2))) # As expected, they're closer to each other
(dist_xz_1 <- sqrt(sum((x - z)^2)))
(dist_yz_1 <- sqrt(sum((y - z)^2)))

# Distance 2: matrix algebra
(dist_xy_2 <- sqrt(crossprod(x - y))) # Point product
(dist_xz_2 <- sqrt(crossprod(x - z)))
(dist_yz_2 <- sqrt(crossprod(y - z)))

# Now to compute all the distances at once (`dist` function)
# Distance 3: Distance between samples (`dist` function)
d <- dist(t(e)) # It's necessary to transpose the full dataset (columns)
class(d)
# Indexes we need to coerce it into a matrix
(dist_xy_3 <- as.matrix(d)[1,2])
(dist_xz_3 <- as.matrix(d)[1,87])
(dist_yz_3 <- as.matrix(d)[2,87])

# """
# It is important to remember that if we run dist on e, 
# it will compute all pairwise distances between genes. 
# This will try to create a 22215×22215 matrix that may 
# crash your R sessions.
# """

# -------------------- Mini test ------------------------
# How many biological replicates are there for hippocampus?
table(tissue)['hippocampus']
# What is the distance between samples 3 and 45?
as.matrix(d)[3, 45]
# What is the distance between gene 210486_at and 200805_at?
x_gene <- e['210486_at',]
y_gene <- e['200805_at',]
sqrt(crossprod(x_gene - y_gene))
sqrt(sum((x_gene - y_gene)^2))
# If I run the command (don't run it!): d = as.matrix(dist(e))
# How many cells (number of rows times number of columns) would this matrix have?
(dim(e)[1]^2) # Without traspose dataset (rows)
# Compute the distance between all pairs of samples: d = dist(t(e))
# How many distances are stored in d? (Hint: What is the length of d)?
floor((dim(e)[2] * (dim(e)[2] - 1))/2)
length(d)
# ----- Because R takes advantage of symmetry: only the lower triangular --------- 
# ----- matrix is stored, thus there are only ncol(e)*(ncol(e)-1)/2 values. ------
