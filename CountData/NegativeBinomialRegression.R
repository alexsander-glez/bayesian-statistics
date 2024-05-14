###################### Negative binomial regression ######################

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

# Import `rjags` library
library("rjags")

##### Model
nb_model_string <- " model{
  ## Likelihood model
  for (i in 1:length(numvisit)){
    numvisit[i] ~ dnegbin(p[i], r)
    p[i] = r / (r + lam[i])
    log(lam[i]) = int + b_badh * badh[i] + b_age * age[i] + b_intx * badh[i] * age[i]
  }

  ## Prior estimations
  int ~ dnorm(0.0, 1.0/1e6) # Non-informative prior
  b_badh ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_age ~ dnorm(0.0, 1.0/1e4) # Non-informative prior
  b_intx ~ dnorm(0.0, 1.0/1e4) # Non-informativa prior
  r ~ dunif(0, 50) # Non-informative prior: r > 0
} "

# To achieve better convergence, let's lengthen the initial adaptation
data_jags <- as.list(badhealth)
nb_params <- c("int", "b_badh", "b_age", "b_intx", "r")

nb_mod <- jags.model(
  textConnection(nb_model_string),
  data = data_jags,
  n.chains = 3,
  n.adapt = 5e3
)
update(nb_mod, 1e3) # burn-in

# Simulated samples from `negative binomial model`
nb_mod_sim <- coda.samples(
  model = nb_mod,
  variable.names = nb_params,
  n.iter = 1e5,
  thin = 5
)

nb_mod_csim <- as.mcmc(do.call(rbind, nb_mod_sim))

### Convergence diagnostics
plot(nb_mod_sim, ask = TRUE)
gelman.diag(nb_mod_sim) # Convergence reached
gelman.plot(nb_mod_sim, ask = TRUE)
autocorr.diag(nb_mod_sim)
autocorr.plot(nb_mod_sim, ask = TRUE)
effectiveSize(nb_mod_sim)
#  b_age    b_badh    b_intx       int         r 
# 6086.652  4536.881  4409.105  6131.043 59251.121
effectiveSize(nb_mod_sim) / 1e5 * 100 # higher dependency between samples
#   b_age  b_badh  b_intx   int     r 
#  6.086%  4.536%  4.409%  6.131% 59.251%
heidel.diag(nb_mod_sim)
raftery.diag(nb_mod_sim)

# Compute DIC
nb_dic = dic.samples(nb_mod, n.iter = 5e3)
# Mean deviance:  4471 
# penalty 4.947 
# Penalized deviance: 4476 

### Model checking
# X: shape(dim) = (1127, 2)
X <- as.matrix(badhealth[,-1]) # column.names = [badh, age]
# X: shape(dim) = (1127, 3)
# colnames(X)[X.columns] = ("badh", "age", "")
X <- cbind(X, with(badhealth, badh * age)) # column.names = [badh, age, -]
# Get the median from each column
(nb_pmed_coef <- apply(nb_mod_csim, 2, median))

# Predicted values
nb_llam_hat <- nb_pmed_coef["int"] + X %*% nb_pmed_coef[c("b_badh", "b_age", "b_intx")] # shape = (1127, 1)
nb_lam_hat <- exp(nb_llam_hat) # shape = (1127, 1)
hist(nb_lam_hat)

### Residual analysis
nb_resid <- badhealth$numvisit - nb_lam_hat
plot(nb_resid) # points are ordered
# plot predicted values vs targets
plot(nb_lam_hat, badhealth$numvisit)
abline(0.0, 1.0)

plot(nb_lam_hat[which(badhealth$badh == 0)],
     nb_resid[which(badhealth$badh == 0)], xlim = c(0, 8),
     ylab = "residuals", xlab = expression(hat(lambda)), ylim = range(nb_resid))
points(nb_lam_hat[which(badhealth$badh == 1)],
       nb_resid[which(badhealth$badh == 1)], col = "red")

# Get residual variance of `badh == 0` and `badh == 1`
var(nb_resid[which(badhealth$badh == 0)]) # var = 7.022656
var(nb_resid[which(badhealth$badh == 1)]) # var = 41.19624

summary(nb_mod_sim)

### Predictive distributions
# What is the posterior probability that the individual with poor health will have more doctor visits?
head(nb_mod_csim) # columns = ["b_age", "b_badh", "b_intx", "int", "r"]
x1 <- c(0, 25, 0) # healthy for Person 1 (good health)
x2 <- c(1, 25, 25) # unhealthy for Person 2 (bad health)
# First, we'll compute the linear part of the predictor
loglam1 <- nb_mod_csim[,"int"] + nb_mod_csim[,c(2, 1, 3)] %*% x1
loglam2 <- nb_mod_csim[,"int"] + nb_mod_csim[,c(2, 1, 3)] %*% x2

# Next we'll apply the inverse link
lam1 <- exp(loglam1)
lam2 <- exp(loglam2)
r1 <- nb_mod_csim[, "r"] # overdispersion parameter
r2 <- nb_mod_csim[, "r"] # overdispersion parameter

# Probabilities
p1 <- r1 / (r1 + lam1)
p2 <- r2 / (r2 + lam2)

# The final step is to use these samples for the $\lambda$ and $r$ parameters for each individual
# and simulate actual number of doctor visits using the likelihood:
(n_sim <- length(lam1))

y1 <- rnegbin(n = n_sim, mu = lam1, theta = r1)
y2 <- rnegbin(n = n_sim, mu = lam2, theta = r2)

plot(
  table(factor(y1, levels = 0:18)) / n_sim,
  pch = 2, ylab = "posterior prob.", xlab = "visits"
)
points(table(y2 + 0.1) / n_sim, col = "red")
