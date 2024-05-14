#### Lesson 10. Poisson regression
## Data
# For an example of Poisson regression, We'll use the badhealth data set from the `COUNT` package in R
if (!require("COUNT")){
  install.packages("COUNT")
  library("COUNT")
}

library("COUNT")
###   Explore data
data("badhealth")
?badhealth
# numvisit -> number of visits to doctor during 1998
# badh -> 1 = patient claims to be in bad health; 0 = not in bad health
# age -> age of patient: 20-60
# Are there any NaN values present in the dataset
any(is.na(badhealth))

# As usual, let's visualize these data
# It follows a Poisson distribution
hist(badhealth$numvisit, breaks = 20) # breaks == bins

plot(
  jitter(log(numvisit)) ~ jitter(age),
  data = badhealth,
  subset = (badh == 0),
  xlab = "age",
  ylab = "log(visits)"
)
points(
  jitter(log(numvisit)) ~ jitter(age),
  data = badhealth,
  subset = (badh == 1),
  col = "red"
)

# It appears that both age and bad health are related to the number of doctor visits. We should include
# model terms for both variables.
# If we believe the age/visits relationship is different between healthy and non-healthy populations,
# we should also include an iteraction term. We will fit the full model here and leave it to you to 
# compare it with the simpler additive model.
library("rjags")

mod_string = " model {
  # Likelihood model (poisson regression)
  for (i in 1:length(numvisit)){
    numvisit[i] ~ dpois(lam[i])
    log(lam[i]) = int + b_badh * badh[i] + b_age * age[i] + b_intx * age[i] * badh[i]
  }
  
  # Prior estimates for each `beta`
  int ~ dnorm(0.0, 1.0 / 1.0e6)
  b_badh ~ dnorm(0.0, 1.0 / 1.0e4)
  b_age ~ dnorm(0.0, 1.0 / 1.0e4)
  b_intx ~ dnorm(0.0, 1.0 / 1.0e4)
} "

set.seed(102)
data_jags <- as.list(badhealth)
params <- c("int", "b_badh", "b_age", "b_intx")
mod <- jags.model(
  textConnection(mod_string),
  data = data_jags,
  n.chains = 3
)
update(mod, 1.0e3) # burn-in
# Model simulation
mod_sim <- coda.samples(
  model = mod,
  variable.names = params,
  n.iter = 5.0e3
)
# Combine all chain simulations in an object
mod_csim <- as.mcmc(do.call(rbind, mod_sim))


### Covergence diagnostics
plot(mod_sim, ask = T)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim) # If exists some autocorr try to increment the number of iterations
effectiveSize(mod_sim)

### compute DIC
dic = dic.samples(mod, n.iter = 1e3)

#### Model checking
# To get a general idea of the model's performance, we can look at predicted values and residuals 
# as usual. Don't forget that we must apply the inverse of the link function to get predictions for 
# $\lambda$
X <- as.matrix(badhealth[,-1])
X <- cbind(X, with(badhealth, badh * age))
head(X)

# Get median coeffiecients for each parameter
(pmed_coef = apply(mod_csim, 2, median))

# log lambda prediction
llam_hat = pmed_coef[1] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)
hist(lam_hat)

# resid <- target - preds 
resid <- badhealth$numvisit - lam_hat
plot(resid) # the data were ordered