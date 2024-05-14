################## POISSON REGRESSION ####################

# badhealth data set from the COUNT package in R
library("COUNT")
data("badhealth")
?badhealth
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

#### Model
# It appears that both age and bad health are related to the number of doctor visits.
# We should include model terms for both variables.
# If we believe the age/visits relationship is different healthy and non-healthy populations, we should
# also include an interaction term. We will fit the full model here and leave it to you to compare it 
# with the simpler additive model.
library("rjags")

mod_string = " model {
  for (i in 1:length(numvisit)){
    numvisit[i] ~ dpois(lam[i])
    # EXPRESSING AS LOGLINK
    log(lam[i]) = int + b_badh * badh[i] + b_age * age[i] + b_intx * badh[i] * age[i]
  }
  
  # Prior estimates
  int ~ dnorm(0.0, 1.0/1e6) # Non-informative prior (high value of variance)
  b_badh ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_age ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_intx ~ dnorm(0.0, 1.0/1e4) # Interaction term, non-informative prior
} "

# Set seed
set.seed(42)
data_jags <- as.list(badhealth) # [numvisit, age, badh]
params <- c("int", "b_badh", "b_age", "b_intx") # Params that need to be calculated
mod <- jags.model(
  textConnection(mod_string), # Model as string format
  data = data_jags, # Variables: numvisit, badh, age
  n.chains = 3 # Maximum number of chains
)
update(mod, 1e3) # burn-in

#### Simulate samples
mod_sim <- coda.samples(
  model = mod,
  variable.names = params, 
  n.iter = 5e3
)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

### Convergence diagnostics
plot(mod_sim, ask = TRUE)
gelman.diag(mod_sim) # Convergence reached -> 1.0
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## Compute DIC
dic = dic.samples(mod, n.iter = 1e3)

### Model checking
# To get a general idea of the model's performance, we can look at predicted values and resisuals.
# Don't forget that we must apply the inverse of the link function to get predictions for \lambda
X = as.matrix(badhealth[,-1]) # All variables, except `numvisit` variable
X = cbind(X, with(badhealth, badh * age)) # cbind -> Column bind (concatenate by columns)
head(X)

# Get median coeffiecients for each parameter
(pmed_coef <- apply(mod_csim, 2, median)) # 2 -> by column
# log lambda prediction
llam_hat <- pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")] # shape = (1127, )
dim(X) # X: shape = (1127, 3)
length(pmed_coef[c("b_badh", "b_age", "b_intx")]) # shape = (3, )
lam_hat <- exp(llam_hat)
hist(lam_hat)

# Analysis of residuals
# resid <- target - preds 
resid <- badhealth$numvisit - lam_hat
plot(resid) # the data were ordered (as exponential behavior)
plot(lam_hat, badhealth$numvisit) # lam_hat: predicted values [vs] badhealth$numvisit: targets 
abline(0.0, 1.0)

# which == where
plot(lam_hat[which(badhealth$badh == 0)], resid[which(badhealth$badh == 0)],
     xlim = c(0, 8), ylab = "residuals", xlab = expression(hat(lambda)), ylim = range(resid))
points(lam_hat[which(badhealth$badh == 1)], resid[which(badhealth$badh == 1)], col = "red")

# It's not surprising that the variability increases for values predicted at higher values since the 
# mean is also the variance in the Poisson distribution. However, observations predicted to have about
# two visits should have variance about two, and observations predicted to have about six visits 
# should have variance about six.
# badh -> In bad health
var(resid[which(badhealth$badh == 0)]) # 7.022435
# badh -> Not in bad health
var(resid[which(badhealth$badh == 1)]) # 41.19662

# Clearly this is not the case with these data. This indicates that either the model fits poorly
# (meaning the covariates don't explain enough of the variablity in the data), or the data are 
# "overdispersed" for the Poisson likelihood we have chosen.
# This is a common issue with the count data. If the data are more variable than the Poisson
# likelihood would suggest, a good alternative is the negative binomial distribution, which we will not pursue here

### Addition model
mod1_string <- " model {
for (i in 1:length(numvisit)){
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh * badh[i] + b_age * age[i]
}
# Prior estimates
int ~ dnorm(0.0, 1.0/1e6) # Non-informative prior
b_badh ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
b_age ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
} "
data1_jags <- as.list(badhealth) # data
params1 <- c("int", "b_badh", "b_age")
# Initialize our model
library("rjags")
mod1 <- jags.model(
  textConnection(mod1_string),
  data = data1_jags,
  n.chains = 3
)
update(mod1, 1e3) # burn-in (stationary state)
mod1_sim <- coda.samples(
  model = mod1,
  variable.names = params1,
  n.iter = 6e3
)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

##### Convergence diagnostic
plot(mod1_sim, ask = TRUE)
gelman.diag(mod1_sim) # Convergence reached
gelman.plot(mod1_sim, ask = TRUE)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim, ask = TRUE)
effectiveSize(mod1_sim)
# compute DIC
dic = dic.samples(mod1, 1e3)

#### Model checking
X <- as.matrix(badhealth[,-1]) # All variables except `numvisit`
head(X)
(pmed1_coef <- apply(mod1_csim, 2, median))
# predicted values
llam1_hat <- pmed1_coef["int"] + X %*% pmed1_coef[c("b_badh", "b_age")]
lam1_hat <- exp(llam1_hat)
# Display predicted values
hist(lam1_hat)
resid1 <- badhealth$numvisit - lam1_hat
plot(resid1) # ths data were ordered
plot(lam1_hat, badhealth$numvisit)
abline(0.0, 1.0)

# 0 -> In bad health
plot(lam1_hat[which(badhealth$badh == 0)], resid1[which(badhealth$badh == 0)],
     xlim = c(0, 8), ylab = "residuals", xlab = expression(hat(lambda)), ylim = range(resid1))
points(lam1_hat[which(badhealth$badh == 1)], resid1[which(badhealth$badh == 1)], col = "red")

var(resid1[which(badhealth$badh == 0)])
var(resid1[which(badhealth$badh == 1)])


