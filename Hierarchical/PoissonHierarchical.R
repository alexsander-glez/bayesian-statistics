########################### Hierarchical Modeling ###########################
##### Data
# Let's fit our hierarchical model for counts of chocolate chips. The data can be found in `cookies.dat`
path_url = "https://raw.githubusercontent.com/007v/Bayesian-Statistics-Techniques-and-Models--University-of-California-Santa-Cruz---Coursera/master/cookies.dat"
dat = read.table(file = path_url, header = TRUE)
head(dat)

# Number of data recorded in each region
table(dat$location)
#  1  2  3  4  5 
# 30 30 30 30 30
# We can also visualize the distribution of chips by location
hist(dat$chips)
boxplot(chips ~ location, data = dat)

##### Prior predictive checks
# We need to select prior distributions for $\alpha$ and $\beta$, the hyperparameters
# governing the gamma distribution for $\lambda$ parameters.
# For location "j", $\lambda_j$ is the expected number of chocolate chips per cookie.
# Hence, $\alpha$ and $\beta$ control the distribution of these means between locations.
# The mean of this gamma distribution will represent the overall mean of number of chips for all cookies.
# The variance of this gamma distribution controls the variability between locations.

n_sim <- 500 # number of simulations
alpha_pri <- rexp(n = n_sim, rate = 1.0/2.0)# draws for $\alpha$ priors ($\lambda_0$ = 1.0/2.0)
beta_pri <- rexp(n = n_sim, rate = 5.0) # draws for $\beta$ priors ($\lambda_0$ = 5.0)
mu_pri <- alpha_pri / beta_pri # prior mean
sig_pri <- sqrt(alpha_pri / beta_pri^2) # standard deviation mean
summary(mu_pri) # Prior mean Summary
#   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0111    3.1529    9.5152   57.0590   33.2173 2212.6771 
summary(sig_pri) # Prior standard deviation summary
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.2295    3.4485    7.7279   40.8186   19.9799 2616.3387 

# After simulating from the priors for $\alpha$ and $\beta$, we can use those
# samples to simulate further down the hierarchy
lam_pri <- rgamma(n = n_sim, shape = alpha_pri, rate = beta_pri)
summary(lam_pri)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    1.167    6.573   54.378   28.934 3834.503 

# Or for a prior predictive reconstruction of the original dataset
(lam_pri <- rgamma(n = 5, shape = alpha_pri[1:5], rate = beta_pri[1:5]))
# lam_pri <- [0.9002759  0.3422614  2.5029962  0.7986137 13.6230912]
(y_pri <- rpois(n = 150, lambda = rep(lam_pri, each = 30)))

# Because these priors have high variance and are somewhat noninformative, they produce unrealistic
# distributions. Still, enough data would overwhelm the prior, resulting in useful posterior distributions.
# Alternatively, we could tweak and simulate from these prior distributions until they adequately represent 
# our prior beliefs.
# Yet another approach would be to re-parametrize the gamma prior, which we'll demonstrate as we fit the model

library("rjags")
head(dat) # Number of chips per cookie by location
pois_mod_sting <- " model{
  # Likelihood model
  for (i in 1:length(chips)){
    # y[i] == chips[i]
    # j = location[i] (1, 2, 3, 4, 5)
    # y[i][j] ~ dpois(lam[location[i]])
    chips[i] ~ dpois(lam[location[i]])
  }
  # Prior estimtations
  for (j in 1:max(location)){
    # max(location) = 5
    lam[j] ~ dgamma(alpha, beta)
  }
  alpha = mu^2 / sig^2 # mu -> mean, sig -> sd
  beta = mu / sig^2 # mu -> mean, sig -> sd
  
  mu ~ dgamma(2.0, 1.0/5.0)
  sig ~ dexp(1.0)
} "
data_jags <- as.list(dat)
pois_params <- c("lam", "mu", "sig")
pois_mod <- jags.model(
  textConnection(pois_mod_sting),
  data = data_jags,
  n.chains = 3
) # initiliaze our model
update(pois_mod, 1e3) # burn-in 
# Draw samples following the hierarchical model
pois_mod_sim <- coda.samples(
  model = pois_mod,
  variable.names = pois_params,
  n.iter = 8e3
)
# Combine all simulations as mcmc 
pois_mod_csim <- as.mcmc(do.call(rbind, pois_mod_sim))

# Compute DIC
pois_dic <- dic.samples(pois_mod, 1e3)

##### Covergence diagnostics
plot(pois_mod_sim, ask = TRUE)
gelman.diag(pois_mod_sim) # Convergence reached
gelman.plot(pois_mod_sim)
autocorr.diag(pois_mod_sim) # Independent samples
autocorr.plot(pois_mod_sim, ask = TRUE)
raftery.diag(pois_mod_sim) # Independence between each variable
heidel.diag(pois_mod_sim)
effectiveSize(pois_mod_sim)
# lam[1]    lam[2]    lam[3]    lam[4]    lam[5]        mu       sig 
# 23229.918 16922.865 22578.710 22541.323 19427.659  9991.534  5662.005

##### Model Checking
# After assessing convergence, we can check the fit via residuals
# With a hierarchical model, there are now two levels of residuals: the observation level 
# and the location mean level.
# To simplify, we'll look at the residuals associated with the posterior means of the parameters

# First, we have OBSERVATION RESIDUALS, based on the estimates of location means.
## observation level residuals
(pm_params <- colMeans(pois_mod_csim))
#  lam[1]    lam[2]    lam[3]    lam[4]    lam[5]        mu       sig 
# 9.277682  6.220389  9.529465  8.946540 11.764175  9.121424  2.097646
# Predicted values
yhat <- rep(pm_params[1:5], each = 30)
pois_resid <- dat$chips - yhat
# Residual plot (it looks ok!)
plot(pois_resid)
# Predicted values vs their resids
plot(jitter(yhat), pois_resid)
var(pois_resid[yhat < 7]) # 6.447126
var(pois_resid[yhat > 11]) # 13.72414

# Also, we can look at how the location means differ from the overall mean $\mu$
# location level residuals
pois_lam_resid <- pm_params[1:5] - pm_params["mu"] # if each mean differ of the global mean
plot(pois_lam_resid)
abline(h = 0, lty = 2)

# We don't see any obvious violations of our model assumptions
summary(pois_mod_sim)


##### Posterior Predictive Simulation
# We can use these posterior samples to get Monte Carlo estimates that interest us from 
# the posterior predictive distribution.
(n_sim <- nrow(pois_mod_csim))
lam_pred <- rgamma(
  n = n_sim,
  shape = pois_mod_csim[, "mu"]^2 / pois_mod_csim[, "sig"]^2,
  rate = pois_mod_csim[, "mu"] / pois_mod_csim[, "sig"]^2
)
hist(lam_pred)
mean(lam_pred) # 9.138441
# Probability that the mean of number of chips per cookie it is more than 15 chips
mean(lam_pred > 15) # 0.01816667

# Using these lambda draws, we can go to the observation level and simulate 
# the number of chips per cookie, which takes into account the uncertainty in lambda
y_pred <- rpois(n = n_sim, lambda = lam_pred) # lam_pred -> lambda with his uncertainty
hist(y_pred)

# Probability that the number of chips per cookie it is more than 15 chips
mean(y_pred > 15) # 0.06120833
hist(dat$chips)

# Finally, we could answer questions like: what is the posterior probability that the next
# cookie produced in Location 1 will have fewer than seven chips?
y_pred1 <- rpois(n = n_sim, lambda = pois_mod_csim[,"lam[1]"])
hist(y_pred1)
mean(y_pred1 < 7) # probability = 0.1904167
