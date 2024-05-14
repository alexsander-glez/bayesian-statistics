# --------------------- Bayesian Linear Regression ------------------
if (!require('car')){
  install.packages('car') # Import the dataset
  library(car)
}
help(Leinhardt) # Some help with Leinhardt dataset
pairs(Leinhardt)
library(ggplot2)
# Omit NA Values
is.na(Leinhardt)
Leinhardt <- na.omit(Leinhardt)
str(Leinhardt)

# Visualizate the data
par(mfrow = c(2,1))
plot(income ~ infant,data = Leinhardt)
plot(infant ~ income, data = Leinhardt)
ggplot(data = Leinhardt) +
  geom_point(mapping = aes(x = income,y = infant),color = 'red',shape = 19)
plot.new()
ggplot(data = Leinhardt) +
  geom_point(mapping = aes(x = infant,y = income),color = 'red',shape = 19)
Leinhardt$region
plot_histogram(as.numeric(Leinhardt$region))

# I we work the data with log-scale
df <- Leinhardt
df$log.income <- log(df$income)
df$log.infant <- log(df$infant)
ggplot(data = df) + 
  geom_point(mapping = aes(x = log.income,y = log.infant))

# 1st Linear Regression (with non-Informative priors)
model.lm <- lm(log.infant ~ log.income,data = df)
summary(model.lm)
plot(resid(model.lm)) # It looks well
plot(predict(model.lm),resid(model.lm)) # It looks like good

# Model using rjags algorithm
model1_string <- " model {
  # Likelihood model
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[i],prec)
    mu[i] = b[1] + b[2]*log.income[i]
  }
  # Prior distributions for beta coefficients
  for (j in 1:2){
    b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative-prior
  }
  # Prior distribution for sigma (standard deviation)
  prec ~ dgamma(5.0/2.0,5.0*9.0/2.0) # sd = 3.0
  sig = sqrt(1/prec) # Squareed root of the variance
} "

# We introduce the data to the model
data1_jags <- list(y = df$log.infant,log.income = df$log.income)
params1 <- c('b','sig')
inits1 <- function(){
  inits = list('b' = rnorm(2,0.0,100),'prec' = rgamma(1,1.0,1.0))
}

# Initiliaze the model
mod1 <- jags.model(textConnection(model1_string),
                   data = data1_jags, inits = inits1,n.chains = 3)
update(mod1,n.iter = 1000) # burn-in

mod1_sim <- coda.samples(model = mod1,
                         variable.names = params1,
                         n.iter = 5.0e3)
mod1_csim <- do.call(rbind,mod1_sim) # Combine multiple chains

# Chck the convergence of the model
plot(mod1_sim)
gelman.plot(mod1_sim)
gelman.diag(mod1_sim) # It converges very well
autocorr.plot(mod1_sim) # b[1] and b[2] are autocorrelated
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

raftery.diag(mod1_sim,q = 0.999,r = 0.001,s = 0.99)

# Predicted values
pm1_params <- colMeans(mod1_csim)
X <- cbind(rep(1.0,length(df$log.income)),data1_jags$log.income)
y1_hat <- drop(X %*% pm1_params[c(1,2)])
resid1 <- df$log.infant - y1_hat
plot(resid1) # It looks like pretty good
abline(h = 0)
plot(y1_hat,resid1); abline(h = 0) # It looks like good

rownames(df)[order(resid1,decreasing = TRUE)[1:5]] # More deviation
summary(mod1_sim)
mean(mod1_csim[,3] > 1.0)
df$oil[rownames(df) == 'Libya']
dic1 <- dic.samples(model = mod1,n.iter = 1.0e3)

# ---------- Linear regression model 2 -------------------
mod2_string <- " model {
  # Likelihood model
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[i],tau)
    mu[i] = b[1] + b[2]*log.income[i] + b[3]*is.oil[i]
  }
  # Response variable {y}
  # Predicted variable {log.income,is.oil}
  ### Prior distribution
  for(j in 1:3){
    b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
  }
  tau ~ dgamma(5.0/2.0,5.0*9.0/2.0)
  sig = sqrt(1.0/tau)
} "

# Incorporing the data
data2_jags <- list(y = df$log.infant,log.income = df$log.income,
                   is.oil = as.numeric(df$oil == 'yes'))
# Model parameters
params2 <- c('b','sig')
 
# Initialization
inits2 <- function(){
  inits = list('b' = rnorm(3,0.0,100),'tau' = rgamma(1,1.0,1.0))
}

# Run the model
mod2 <- jags.model(textConnection(mod2_string),
                   data = data2_jags,
                   inits = inits2,n.chains = 3)
update(mod2,1.0e3) # burn-in
mod2_sim <- coda.samples(model = mod2,
                         variable.names = params2,
                         n.iter = 7.0e3)
mod2_csim <- do.call(rbind,mod2_sim) # Combine multiple chains

# Check the convergence for the model 2
plot(mod2_sim)
gelman.diag(mod2_sim) # It converges well
gelman.plot(mod2_sim)
autocorr.plot(mod2_sim)
autocorr.diag(mod2_sim)
effectiveSize(mod2_sim)
raftery.diag(mod2_sim)

# Prediction (model 2) [Additional covariates]
(pm2_params <- colMeans(mod2_csim))
X2 <- cbind(rep(1.0,nrow(df)),df$log.income,as.numeric(df$oil == 'yes'))
yhat2 <- drop(X2 %*% pm2_params[1:3])
resid2 <- df$log.infant - yhat2  
plot(resid2); abline(h = 0,col = 'red') # It looks like good
plot(yhat2,resid2) # It looks like good too

dic2 <- dic.samples(model = mod2,n.iter = 1.0e3)

# ---------- Linear regression model 3 ------------------
curve(dnorm(x,mean = 0.0,sd = 1.0),lty = 1,xlim = c(-4,4))
curve(dt(x,df = 1),col = 'red',type = 'l',add = TRUE)

mod3_string <- " model {
  # Likelihood model
  for(i in 1:length(y)){
    y[i] ~ dt(mu[i],tau,df)
    mu[i] = b[1] + b[2]*log.income[i] + b[3]*is.oil[i]
  }
  # Prior estimation
  for(j in 1:3){
    b[j] ~ dnorm(0.0,1.0/1.0e6) # Non-informative prior
  }
  df = nu + 2.0 # For it doesn't exist indeterminations
  nu ~ dexp(1.0)
  tau ~ dgamma(5.0/2.0,5.0*9.0/2.0)
  sig = sqrt(1/tau*df/(df - 2))
} "

# We incorporate the data
data3_jags <- list(y = df$log.infant,log.income = df$log.income,
                   is.oil = as.numeric(df$oil == 'yes'))
# Variable names
params3 <- c('b','sig')

# Inits
inits3 <- function(){
  inits = list('b' = rnorm(3,0.0,100.0),'tau' = rgamma(1,1.0,1.0))
}

# Run the model
mod3 <- jags.model(textConnection(mod3_string),
                   data = data3_jags,inits = inits3,n.chains = 3)
update(mod3,1.0e3) # burn-in
mod3_sim <- coda.samples(model = mod3,
                         variable.names = params3,
                         n.iter = 5.0e3)
mod3_csim <- do.call(rbind,mod3_sim) # Combine multiple chains

# Check the convergence (model 3)
plot(mod3_sim)
gelman.diag(mod3_sim) # It converges well
gelman.plot(mod3_sim)
autocorr.plot(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)

# Predictive model
(pm3_params <- colMeans(mod3_csim))
(X3 <- X2); y <- df$log.infant;
y3_hat <- drop(X3 %*% pm3_params[1:3]) 
resid3 <- y - y3_hat
plot(resid3); abline(h = 0,col = 'red')
plot(y3_hat,resid3); abline(h = 0,col = 'red')
dic3 <- dic.samples(mod3,n.iter = 1.0e3)

dic1
dic2 # The better model
dic3
