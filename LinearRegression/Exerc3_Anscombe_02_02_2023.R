if(require('car')){
  install.packages('car')
  library(car)
  data("Anscombe")
}
help("Anscombe")
str(Anscombe); df <- Anscombe
pairs(df)

# Adjust the first model
model.lm <- lm(education ~ income + young + urban,data = df)
names(df)
summary(model.lm) # Reference (non-informative) Bayesian linear model
coef(model.lm)

## Initial guess of variance based on overall
## variance of education variable. Uses low prior
## effective sample size. Technically, this is not
## a true 'prior', but it is not very informative.

model_string <- " model {
  # Likelihood model
  for(i in 1:length(education)){
    education[i] ~ dnorm(mu[i],tau)
    mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i] 
  }
  # Prior distribution for beta coefficients
  b0 ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
  for(j in 1:3){
    b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
  }
  tau ~ dgamma(1/2.0,1*1500.0/2.0)
  sig2 = 1/tau
} "
data_jags <- as.list(df)
params <- c('b','sig2')
inits <- function(){
  inits <- list('b0' = rnorm(1,0.0,100.0),'b' = rnorm(3,0.0,100.0),
                'tau' = rgamma(1,1.0,1.0))
}
# Run the model
mod <- jags.model(textConnection(model_string),
                  data = data_jags,inits = inits,n.chains = 4)
update(mod,1.0e3) # burn-in
mod_sim <- coda.samples(model = mod,
                        variable.names = params,
                        n.iter = 100.0e3)
mod_csim <- do.call(rbind,mod_sim) # Combine multiple chains
mean(mod_csim[,1] > 0.0)
# ----------------- MCMC Convergence ----------------
plot(mod_sim)
gelman.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
(dic1 <- dic.samples(mod,n.iter = 1.0e5))
plot(resid(model.lm))
plot(model.lm)
summary(mod_sim)

# ---------------------- Model 2----------------------------------
plot_correlation(df)

model2_string <- " model {
    # Likelihood model
    for(i in 1:length(education)){
      education[i] ~ dnorm(mu[i],tau)
      mu[i] = b0 + b[1]*income[i] + b[2]*young[i] 
    }
    # Prior distribution for beta coefficients
    b0 ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
    for(j in 1:2){
      b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
    }
    tau ~ dgamma(1/2.0,1*1500.0/2.0)
    sig2 = 1/tau
  } "

data2_jags <- list(education = df$education,income = df$income,
                      young = df$young)
params2 <- c('b','sig2')
inits2 <- function(){
  inits <- list('b0' = rnorm(1,0.0,100.0),'b' = rnorm(2,0.0,100.0),
                'tau' = rgamma(1,1.0,1.0))
}
# Run the model
mod2 <- jags.model(textConnection(model2_string),
                  data = data2_jags,inits = inits2,n.chains = 4)
update(mod2,1.0e4) # burn-in
mod2_sim <- coda.samples(model = mod2,
                        variable.names = params2,
                        n.iter = 100.0e3)
mod2_csim <- do.call(rbind,mod2_sim) # Combine multiple chains

# ----------------- MCMC Convergence ----------------
plot(mod2_sim)
gelman.diag(mod2_sim)
effectiveSize(mod2_sim)
(dic2 <- dic.samples(mod2,n.iter = 1.0e5))
plot(resid(model.lm))
plot(model.lm)


# ------------------------- Model 3 -------------------------------
plot_correlation(df)
model3_string <- " model {
    # Likelihood model
    for(i in 1:length(education)){
      education[i] ~ dnorm(mu[i],tau)
      mu[i] = b0 + b[1]*income[i] + b[2]*income[i]*young[i] 
    }
    # Prior distribution for beta coefficients
    b0 ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
    for(j in 1:2){
      b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
    }
    tau ~ dgamma(1/2.0,1*1500.0/2.0)
    sig2 = 1/tau
  } "

data3_jags <- list(education = df$education,income = df$income,
                   young = df$young)
params3 <- c('b','sig2')
inits3 <- function(){
  inits <- list('b0' = rnorm(1,0.0,100.0),'b' = rnorm(2,0.0,100.0),
                'tau' = rgamma(1,1.0,1.0))
}
# Run the model
mod3 <- jags.model(textConnection(model3_string),
                   data = data3_jags,inits = inits3,n.chains = 4)
update(mod3,1.0e4) # burn-in
mod3_sim <- coda.samples(model = mod3,
                         variable.names = params3,
                         n.iter = 100.0e3)
mod3_csim <- do.call(rbind,mod3_sim) # Combine multiple chains

# ----------------- MCMC Convergence ----------------
plot(mod3_sim)
gelman.diag(mod2_sim)
effectiveSize(mod2_sim)
(dic3 <- dic.samples(mod3,n.iter = 1.0e5))
plot(resid(model.lm))
plot(model.lm)