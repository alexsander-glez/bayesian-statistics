# Prior t-distribution
curve(dt(x,df = 1),col = 'blue',type = 'l',xlim = c(-3.0,3.0))

# Target distribution (Normal likelihood (known variance) and t-distribution prior)
# log-scale -> Numerically stable (likelihoods product)
lg <- function(n,ybar,mu){
  mu2 <- mu^2;
  return(n*(mu*ybar - mu2/2.0) - log(1 + mu2))
}

# Data recollection (Percent change in total personnel from last year to this year)
y <- c(1.2,1.4,-0.5,0.3,0.9,2.3,1.0,0.1,1.3,1.9); n <- length(y)
ybar <- mean(y);
hist(y,freq = FALSE,xlim = c(-1.0,3.0)) # Histogram of the data
curve(dt(x,df = 1),col = 'blue',type = 'l',add = TRUE) # Prior for mu
points(y,rep(0,n),pch = 1) # Individual data points

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


post1 <- mh(n = n,ybar = ybar,n_iter = 1.0e3,mu_init = 0,cand_sd = 3.0)
post1$accpt # Accept ratio below 23% 
library('coda')
traceplot(as.mcmc(post1$mu))
list(post1)

post2 <- mh(n = n,ybar = ybar,n_iter = 1.0e3,mu_init = 0,cand_sd = 0.05)
post2$accpt # Accept rate above 50%
traceplot(as.mcmc(post2$mu)) # It doesn't cover the space well
list(post2) 

post3 <- mh(n = n,ybar = ybar,n_iter = 1.0e3,mu_init = 0,cand_sd = 0.9)
post3$accpt # 23% < Accept rate < 50% (Good)
traceplot(as.mcmc(post3$mu))
list(post3)
# Let's see what happens if we initialize the chain at some far-off value
post4 <- mh(n = n,ybar = ybar,n_iter = 1.0e3,mu_init = 30.0,cand_sd = 0.9)
traceplot(as.mcmc(post4$mu))
post4$accpt

# Burn-in
post4$mu_keep <- post4$mu[-c(1:100)] # Discard the first 100 samples
traceplot(as.mcmc(post4$mu_keep))
plot(density(post4$mu_keep,adjust = 2.0),main = '',xlim = c(-1.0,3.0),
     xlab = expression(mu)) # Density estimate of the posterior
curve(dt(x,df = 1),lty = 2,col = 'blue',add = TRUE) # Prior for mu
points(ybar,0,pch = 19) # Sample mean
points(y,rep(0,n),pch = 1)

# Approximation to the true posterior in red
curve(0.017*exp(lg(mu = x,n = n,ybar = ybar)),from = -1.0,to = 3.0,add = TRUE,col = 'red')

# ------------------- JAGS Sampler -------------------
install.packages('rjags')
library(rjags)
# Data recollection (Percent change in total personnel from last year to this year)
y <- c(1.2,1.4,-0.5,0.3,0.9,2.3,1.0,0.1,1.3,1.9); n <- length(y)

# 1) Specify the model
mod_string <- " model {
  # Likelihood model
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu,1.0/sig2)
  }
  # Prior distributions
  mu ~ dt(0,1.0/1.0,1) # Location, inverse scale, degrees of freedom
  sig2 = 1.0
}"

# 2) Set up the model
data_jags <- list(y = y)
params <- c('mu')
inits <- function(){
  list('mu' = 0)
} # Optional (and fixed)
mod <- jags.model(textConnection(mod_string),data = data_jags,inits = inits)

# 3) Run the MCMC Sampler
update(mod,n.iter = 1.0e3) # burn-in
mod_sim <- coda.samples(model = mod,
                        variable.names = params,
                        n.iter = 5.0e3)

# 4) Post Processing
library('coda')
plot(mod_sim)
summary(mod_sim)
