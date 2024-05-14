####################### Poisson regression #######################

# badhealth data set from the COUNT package in R
library("COUNT")
data("badhealth") # `badhealth` dataset
head(badhealth)

# Are there any NAN values into dataset?
any(is.na(badhealth)) # -> FALSE
# Let's visualize the data (COUNT -> histogram)
hist(badhealth$numvisit, breaks = 20) # breaks == bins
# Relation between age and number of visits
plot(jitter(log(numvisit)) ~ jitter(age), data = badhealth, 
     subset = badh == 0, xlab = "age", ylab = "log(visits)") # In bad health
points(jitter(log(numvisit)) ~ jitter(age), data = badhealth,
       subset = (badh == 1), col = "red") # Not in bad health

library("rjags")

##### Model
pois_model_string <- " model {
  ## Likelihood model
  for (i in 1:length(numvisit)){
    numvisit[i] ~ dpois(lam[i])
    log(lam[i]) = int + b_badh * badh[i] + b_age * age[i] + b_intx * badh[i] * age[i]
  }
  ## Prior estimations
  int ~ dnorm(0.0, 1.0/1e6) # Non-informative prior
  b_badh ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_age ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_intx ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
} "
data_jags <- as.list(badhealth)
pois_params <- c("int", "b_badh", "b_age", "b_intx")

pois_mod <- jags.model(
  textConnection(pois_model_string),
  data = data_jags,
  n.chains = 3
) # compile model
update(pois_mod, 1e3) # burn-in

# Time to update our model, and generate posteriors for our two regression parameters
pois_mod_sim <- coda.samples(
  model = pois_mod,
  variable.names = pois_params,
  n.iter = 1e5
) # increase chain length from default

# Now we can check that JAGS gets the mean of the posteriors of our parameters roughly right (close to 1)
# and compare them to the ML estimates we get from `glm`:
pois_mod_csim <- as.mcmc(do.call(rbind, pois_mod_sim))
summary(pois_mod_sim)

### Poisson convergence diagnostics
plot(pois_mod_sim, ask = TRUE)
gelman.diag(pois_mod_sim) # Convergence reached
gelman.plot(pois_mod_sim, ask = TRUE)
autocorr.diag(pois_mod_sim)
autocorr.plot(pois_mod_sim, ask = TRUE)
effectiveSize(pois_mod_sim)
#   b_age   b_badh   b_intx      int 
# 4592.557 3823.360 3647.577 4682.585
# We should also assess the posteriors' stationarity (Heidelberg-Welch convergence diagnostic)
heidel.diag(pois_mod_sim) # Seems ok. Let's also check that our chain's length is satisfactory
raftery.diag(pois_mod_sim)
# This corresponds to only (4.81%, 3.82%, 3.74%, 5.02%) of our draws, signaling strong dependency 
# in our chain.

# Compute DIC
p_dic <- dic.samples(pois_mod, n.iter = 5e3)
# Mean deviance:  5630 
# penalty 3.975 
# Penalized deviance: 5634

### Model checking
X <- as.matrix(badhealth[,-1]) # column.names = [badh, age]
# dim(X): shape = (1127, 3)
X <- cbind(X, with(badhealth, badh * age)) # column.names = [badh, age, -]

# Get median of each column from simulated samples
(pois_pmed_coef <- apply(pois_mod_csim, 2, median))
#    b_age       b_badh       b_intx          int 
# 0.008501159  1.566098030 -0.010860321  0.346614901

# predicted values
# Get median coeffiecients for each parameter
pois_llam_hat <- pois_pmed_coef["int"] + X %*% pois_pmed_coef[c("b_badh", "b_age", "b_intx")] # shape = (1127, 1)
pois_lam_hat <- exp(pois_llam_hat) # shape = (1127, 1)
hist(pois_lam_hat)

### Residual analysis
pois_resid <- badhealth$numvisit - pois_lam_hat
plot(pois_resid) # data is ordered
# Plot predicted vs targets
plot(pois_lam_hat, badhealth$numvisit)
abline(0.0, 1.0)

# Not in bad health (~sick)
plot(pois_lam_hat[which(badhealth$badh == 0)], 
     pois_resid[which(badhealth$badh == 0)], xlim = c(0, 8), ylab = "residuals",
     xlab = expression(hat(lambda)), ylim = range(pois_resid))
points(pois_lam_hat[which(badhealth$badh == 1)],
       pois_resid[which(badhealth$badh == 1)], col = "red") # In bad health (sick)

# It is not surprising that the variability increases for values predicted at higher values
# since the mean is also the variance in the Poisson distribution.
var(pois_resid[which(badhealth$badh == 0)]) # var = 7.022528
var(pois_resid[which(badhealth$badh == 1)]) # var = 41.19612

# summary(pois_mod_sim)

### Predictive distributions
# What is the posterior probability that the individual with poor health will have more doctor visits?
head(pois_mod_csim) # columns = ["b_age", "b_badh", "b_intx", "int"]
x1 <- c(0, 25, 0) # healthy for Person 1 (good health)
x2 <- c(1, 25, 25) # unhealthy for Person 2 (bad health)
# First, we'll compute the linear part of the predictor
loglam1 <- pois_mod_csim[,"int"] + pois_mod_csim[,c(2, 1, 3)] %*% x1
loglam2 <- pois_mod_csim[,"int"] + pois_mod_csim[,c(2, 1, 3)] %*% x2

# Next we'll apply the inverse link
lam1 <- exp(loglam1)
lam2 <- exp(loglam2)

# The final step is to use these samples for the $\lambda$ parameter for each individual
# and simulate actual number of doctor visits using the likelihood:
(n_sim <- length(lam1))

y1 <- rpois(n = n_sim, lambda = lam1)
y2 <- rpois(n = n_sim, lambda = lam2)

plot(
  table(factor(y1, levels = 0:18)) / n_sim,
  pch = 2, ylab = "posterior prob.", xlab = "visits"
)
points(table(y2 + 0.1) / n_sim, col = "red")
