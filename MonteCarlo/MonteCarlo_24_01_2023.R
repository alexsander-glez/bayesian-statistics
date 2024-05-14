# $\theta$ ~ Gamma(a,b) with a = 2 & b = 1/3
set.seed(32) # Initializes the random number generator so we can replicate these results
# To get different numbers, change the random seed
m <- 1.0e4; a <- 2; b <- 1.0/3

# To simulate m independent samples, use the "rgamma" function
theta <- rgamma(n = m,shape = a,rate = b)
# We can plot a histogram of the generated data, and compare that to the true density
par(mfrow = c(1,1))
hist(theta,freq = FALSE)
# Gamma distribution with that hyperparameters
curve(dgamma(x,shape = a,rate = b),col = 'blue',add = TRUE)

# To find our Monte Carlo estimation, let's take the average of our sample
(true.mean <- a/b)
(approx.mean <- sum(theta)/m)

gamma.means <- function(m,a = 2,b = 1.0/3){
  true.mean <- a/b
  theta <- rgamma(n = m,shape = a,rate = b)
  approx.mean <- mean(theta)
  return(c(true.mean,approx.mean))
}
gamma.means(m = 1.0e4)

# How about the variance of theta?
var(theta) # Sample variance
a/b^2 # True variance of gamma

# We can also approximate the probability that theta < 5
mean(theta < 5)
pgamma(q = 5.0,shape = a,rate = b) # The true value of P(theta < 5)
boxplot(theta)
quantile(theta,probs = 0.90)
qgamma(p = 0.90,shape = a,rate = b) # The true 90% quantile for the data

# ------------------------- Monte Carlo Error -------------------------
# We use the CLT to approximate how accurate our Monte Carlo estimates are
# Monte Carlo standard deviation
se <- sqrt(var(theta)/m)
2.0*se # We are reasonably confident that the Monte Carlo estimate is no more than this 
# far from the truth
c(approx.mean - 2.0*se,approx.mean + 2.0*se) # 95% Confident (two standard deviations)

# These numbers give us a reasonable range for the quantity we are estimating with Monte Carlo
# The same applies for other Monte Carlo estimates, like probability that theta < 5
se <- sqrt(var(theta < 5)/m)
c(approx.mean - 2.0*se,approx.mean + 2.0*se)
2.0*se # We are reasonably confident that the Monte Carlo estimate is no more than this
# far from the truth

# --------------------- Marginalization ---------------------
# Let's also do the second example of simulation a hierarchical model
# y|phi ~{iid}~ Bin(10,phi)
# phi ~ Beta(2.0,2.0)
# Simulate from this joint distribution
m <- 1.0e4; alpha <- 2.0; beta <- 2.0;

y <- numeric(m)
phi <- numeric(m) # Create the vectors we fill in with simulations

phi <- rbeta(n = m,shape1 = alpha,shape2 = beta)
y <- rbinom(n = m,size = 10,prob = phi)

# If we are interested only in the marginal distribution of y
mean(y)
plot(prop.table(table(y)),ylab = "P(y)", main = "Marginal distribution of y")
head(y)
table(y) # Count the elements and reagrupe in one group
sum(prop.table(table(y))) # Sum of all probabilities from each element
