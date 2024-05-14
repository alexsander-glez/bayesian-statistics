# ------------------------- GIBBS SAMPLER IN R ----------------------
# Data -> The percent change of personnel from last year to this year for 10 companies
# We'll still use a normal likelihood, but now we'll relax the assumption that we
# know the variance of growth between companies (estimate the variance)

# Prior distribution estimate (before observe the data)
prior <- list()
prior$mu_0 <- 0.0
prior$sig2_0 <- 1.0
curve(dnorm(x,mean = prior$mu_0,sd = sqrt(prior$s2_0)),lty = 2,xlim = c(-1.0,3.0))

# Collected data (for 10 companies)
# Explore the data
y <- c(1.2,1.4,-0.5,0.3,0.9,2.3,1.0,0.1,1.3,1.9)
n <- length(y); ybar <- mean(y)
hist(y,freq = FALSE,xlim = c(-1.0,3.0)) # Histogram of the data
points(y,rep(0,n),pch = 1) # Data points
points(ybar,0,pch = 19) # Sample mean
curve(dnorm(x,mean = prior$mu_0,sd = sqrt(prior$s2_0)),lty = 2,add = TRUE)

# Posterior distribution
# Full conditional distribution for mu (update mu)
update_mu <- function(n,ybar,sig2,mu.0,s0.2){
  s1.2 <- 1.0/(n/sig2 + 1.0/s0.2)
  mu.1 <- s1.2*(n*ybar/sig2 + mu.0/s0.2)
  rnorm(n = 1,mean = mu.1,sd = sqrt(s1.2))
}

# Full conditional distribution for sig2 (update sig2)
update_sig2 <- function(n,y,mu,nu.0,beta.0){
  nu.1 <- nu.0 + n/2.0
  sumsq <- sum((y - mu)^2)
  beta.1 <- beta.0 + sumsq/2.0
  out_gamma <- rgamma(n = 1,shape = nu.1,rate = beta.1)
  1.0/out_gamma # Reciprocal of gamma random variable is distributed inv-gamma
}

# --------------------- Gibbs sampling -----------------
gibbs <- function(y,n_iter,init,prior){
  ybar <- mean(y)
  n <- length(y)
  
  # 1) Initialize
  mu_out <- numeric(n_iter)
  sig2_out <- numeric(n_iter)
  
  mu_now <- init$mu_0
  
  # 2) for i = 1,...,m, repeat:
  for (i in 1:n_iter){
    sig2_now <- update_sig2(n = n,y = y,mu = mu_now,nu.0 = prior$nu_0,beta.0 = prior$beta_0)
    mu_now <- update_mu(n = n,ybar = ybar,sig2 = sig2_now,mu.0 = prior$mu_0,s0.2 = prior$sig2_0)
  
    sig2_out[i] <- sig2_now
    mu_out[i] <- mu_now
  }
  cbind(mu = mu_out,sig2 = sig2_out)
}

# Prior estimates (for sig2)
prior$n_0 <- 2.0 # Effective prior sample size {alpha_0}
prior$s2_0 <- 1.0 # Prior point estimate for sig2
prior$nu_0 <- prior$n_0/2.0 # Prior paramater for inverse-gamma
prior$beta_0 <- prior$n_0*prior$s2_0/2.0 # Prior parameter for inverse-gamma

# Finally, we can initialize and run the sampler
init <- list()
init$mu_0 <- 0.0
post <- gibbs(y = y,n_iter = 1.0e3,init = init,prior = prior)

# Post Processing
library('coda')
plot(as.mcmc(post))
summary(as.mcmc(post))
