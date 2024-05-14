# --------------------- ANOVA ANALYSIS -------------------
library('rjags')
library('coda')
data("PlantGrowth") # Import dataset
df <- PlantGrowth
str(df)
# Boxplot (weight ~ group)
boxplot(weight ~ group,data = df)
help("PlantGrowth")

# Refeerence (noninformative) Bayesian Anova Analysis
lm.mod <- lm(weight ~ group,data = df)
summary(lm.mod)
ctrl.mean.lm <- 5.0320
trt1.mean.lm <- ctrl.mean.lm - 0.3710
trt2.mean.lm <- ctrl.mean.lm + 0.4940
cbind(ctrl.mean.lm,trt1.mean.lm,trt2.mean.lm)
lm.aov <- aov(lm.mod) # 95% confident that at least one mean is different
TukeyHSD(lm.aov,conf.level = 0.99)

# -------------------- RJAGS MODEL ANOVA ----------------------
model_string <- " model {
  # Likelihood model
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[grp[i]],tau)
  }
  # Prior distributions (estimations)
  for(j in 1:3){
    mu[j] ~ dnorm(0.0,1.0/1.0e6) # Noninformative prior {flat prior}
  }
  tau ~ dgamma(3/2.0,3*1.0/2.0)
  sig = sqrt(1/tau)
} "

# Incorporing the data
data_jags <- list(y = df$weight,grp = as.numeric(df$group))
params <- c('mu','sig')
inits <- function(){
  inits <- list('mu' = rnorm(3,0.0,100.0),'tau' = dgamma(1,1.0,1.0))
}
# Run the model
mod <- jags.model(textConnection(model_string),
                  data = data_jags, inits = inits, n.chains = 3)
update(mod,1.0e3)
mod_sim <- coda.samples(model = mod,
                        variable.names = params,
                        n.iter = 6.0e3)
mod_csim <- as.mcmc(do.call(rbind,mod_sim)) # Combine multiple chains

# -------------------- MCMC Convergence ------------------------
plot(mod_sim)
gelman.diag(mod_sim) # It converges really well
gelman.plot(mod_sim)
autocorr.diag(mod_sim) # Is not any autocorrelation between samples
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

# Predictive model
(pm_params <- colMeans(mod_csim))
(y_hat <- pm_params[1:3][data_jags$grp]) 
y <- data_jags$y
resid1 <- y - y_hat
plot(resid1) # It looks pretty good
plot(y_hat,resid1) # Control (not treatment) has most variance

TukeyHSD(lm.aov)
HPDinterval(mod_csim,prob = 0.95)
HPDinterval(mod_csim[,3] - mod_csim[,2],0.95)
mean(mod_csim[,3] > 1.1*mod_csim[,1])
