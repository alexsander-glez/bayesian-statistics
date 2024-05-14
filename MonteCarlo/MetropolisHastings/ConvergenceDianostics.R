# ------------------ CONVERGENCE DIAGNOSTICS ------------------
# Has our simulated Markov chain converged to its stationary distribution yet?
y <- c(1.2,1.4,-0.5,0.3,0.9,2.3,1.0,0.1,1.3,1.9); n <- length(y)
ybar <- mean(y);

# TRACE PLOTS
set.seed(61)
# A traceplot shows the history of a parameter value across iterations of the chain
# (Where the chain has been exploring)
post0 <- mh(n = n,ybar = ybar,n_iter = 10e3,mu_init = 0.0,cand_sd = 0.9)
library('coda')
traceplot(as.mcmc(post0$mu[-c(1:500)])) # If the chain is stationary, it should
# not be showing any long-term trends

set.seed(61)
post1 <- mh(n = n,ybar = ybar,n_iter = 1e3,mu_init = 0.0,cand_sd = 0.04)
traceplot(as.mcmc(post1$mu[-c(1:500)]))
# Is this the case, you need to run the chain many more iterations
set.seed(61)
post2 <- mh(n = n,ybar = ybar,n_iter = 100e3,mu_init = 0.0,cand_sd = 0.04)
traceplot(as.mcmc(post2$mu)) # The chain appear to have converged at this much larger 
# time scale

# MONTE CARLO EFFECTIVE SAMPLE SIZE
# Level of autocorrelation in each
# How linearly depend the current value of the chain is on past values (called lags)
autocorr.plot(as.mcmc(post0$mu))
autocorr.diag(as.mcmc(post0$mu))
autocorr.plot(as.mcmc(post1$mu)) # High autocorrelation
autocorr.diag(as.mcmc(post1$mu))

# Autocorrelation is important because it tells us how much information is available
# in our Markov Chain
str(post2)
(efss1 <- effectiveSize(as.mcmc(post2$mu))) # effecctive sample size (high autocorrelated chain)
(efss2 <- effectiveSize(as.mcmc(post0$mu)))
round(efss2/efss1,0)

# Thin out the samples until correlation is essentially 0.
# This will leave you with approximately independent samples
# The number of samples remaining is similar to the effective sample size
autocorr.plot(as.mcmc(post2$mu),lag.max = 500)
thin_interval <- 400 # how far apart the iterations are for autocorrelation to be essentially 0
thin_indx <- seq(from = thin_interval,to = length(post2$mu),by = thin_interval)
head(thin_indx)

autocorr.plot(as.mcmc(post2$mu[thin_indx]))
effectiveSize(as.mcmc(post2$mu[thin_indx]))
length(thin_indx)

post2mu_thin <- post2$mu[thin_indx]
traceplot(as.mcmc(post2$mu))
effectiveSize(as.mcmc(post2$mu)) # Effective sample size ~ 370
length(post2$mu) # 100,000 iterations
traceplot(as.mcmc(post2mu_thin))
effectiveSize(as.mcmc(post2mu_thin)) # Effective sample size ~ 250
length(post2mu_thin) # 250 iterations

# post2$mu give us more information about stationary distribution
# instead, it has more number of iterations (~ 400 times)
autocorr.plot(as.mcmc(post2mu_thin),lag.max = 10)

str(post0) # it has 10,000 iterations
effectiveSize(as.mcmc(post0$mu)) # ~ 2500 independent MC samples
?effectiveSize

# Raftery-Lewis diagnostic
raftery.diag(as.mcmc(post0$mu)) # We need 3746 independent MC samples
# to calculate 95% posterior interval (reliable estimate)
raftery.diag(post0$mu,q = 0.005,r = 0.001,s = 0.95)
?raftery.diag

## BURN-IN
# We have also seen how the initial value of the chain can affect how quickly the chain
# converges
set.seed(62)
post3 <- mh(n = n,ybar = ybar,n_iter = 500,mu_init = 10.0,cand_sd = 0.3)
traceplot(as.mcmc(post3$mu)) # It takes some time to reach the stationary distribution
traceplot(as.mcmc(post3$mu[-c(1:100)])) # Burn-in the first 100 samples (it reaches the stationary distribution)
# You should always discard early iterations that do not appear te be coming from the stationary distribution

## MULTIPLE CHAINS, GELMAN-RUBIN
# If we want to be more confident that we have converged to the true stationary distribution
# We can simulate multiple chains, each with different starting value
set.seed(61)
nsim <- 500; mu_init <- c(15.0,-5.0,7.0,23.0,-17.0)
cand_sd <- c(0.4,0.4,0.1,0.5,0.4)
post1 <- mh(n = n,ybar = ybar,n_iter = nsim,mu_init = mu_init[1],cand_sd = cand_sd[1])
post2 <- mh(n = n,ybar = ybar,n_iter = nsim,mu_init = mu_init[2],cand_sd = cand_sd[2])
post3 <- mh(n = n,ybar = ybar,n_iter = nsim,mu_init = mu_init[3],cand_sd = cand_sd[3])
post4 <- mh(n = n,ybar = ybar,n_iter = nsim,mu_init = mu_init[4],cand_sd = cand_sd[4])
post5 <- mh(n = n,ybar = ybar,n_iter = nsim,mu_init = mu_init[5],cand_sd = cand_sd[5])

post1$accpt
post2$accpt
post3$accpt
post4$accpt
post5$accpt

pmc <- mcmc.list(as.mcmc(post1$mu),as.mcmc(post2$mu),as.mcmc(post3$mu),
                 as.mcmc(post4$mu),as.mcmc(post5$mu))

str(pmc)
traceplot(pmc) # All chains are exploring the stationary (posterior) distribution
# Gelman and Rubindiagnostics
# Variability within chains, comparing that to the variability between chains
gelman.diag(pmc)
gelman.plot(pmc) # We have probably reached convergence

## MONTE CARLO ESTIMATION
# If we are reasonably confident that our Markov Chain has converged
nburn <- 1000 # Discard early iterations
post0$mu_keep <- post0$mu[-c(1:nburn)]
summary(as.mcmc(post0$mu_keep))
mean(post0$mu_keep > 1.0) # Posterior probability that mu > 1.0

# Auxiliary functions (Metropolis-Hastings algorithm)

mh <- function(n,ybar,n_iter,mu_init,cand_sd){
  # Random Walk Metropolis-Hastings algorithm
  # 1) Select an initial value (mu_init)
  mu_out <- numeric(n_iter)
  accpt <- 0
  mu_now <- mu_init
  lg_now <- lg(n = n,ybar = ybar,mu = mu_now)
  
  # 2) For i = 1,...,m:
  for(i in 1:n_iter){
    # a) Propose condidate (mu*) {to be other state as mu_{i-1}}
    mu_cand <- rnorm(n = 1,mean = mu_now,sd = cand_sd)
    lg_cand <- lg(n = n,ybar = ybar,mu = mu_cand)
    lalpha <- lg_cand - lg_now # log of acceptance ratio
    # b) Calculate alpha ratio
    alph <- exp(lalpha) # Accept or Reject ratio
    
    # c) Accept or reject the proposal candidates
    u <- runif(1)
    if (alph > u){
      mu_now = mu_cand # Accept the draw sample
      accpt <- accpt + 1
      lg_now <- lg_cand
    }
    # Collect results
    mu_out[i] <- mu_now # Reject the proposal draw
  }
  ## Return a list of output
  list(mu = mu_out,accpt = accpt/n_iter)
}

# Target distribution (Normal likelihood (known variance) and t-distribution prior)
# log-scale -> Numerically stable (likelihoods product)
lg <- function(n,ybar,mu){
  mu2 <- mu^2;
  return(n*(mu*ybar - mu2/2.0) - log(1 + mu2))
}

