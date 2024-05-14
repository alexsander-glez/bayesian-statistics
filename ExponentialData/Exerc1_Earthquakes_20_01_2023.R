# Exponential Data (waiting times)
lamb_ <- seq(from = 1, to = 20, by = 1)
lamb <- 1/lamb_

# ------------------ Prior Disribution --------------------
# Prior mean = 1/20; Effective sample size = 1;
alpha <- 1; beta <- 20;
#   Prior Distribution (Gamma Distribution)
plot(lamb,dgamma(lamb,shape = alpha,rate = beta),type = 'l',col = 'black')

# ------------------ Data ---------------------------
y <- c(12,15,8,13.5,25); n <- length(y)

# -------------- Posterior Distribution ------------------
alpha.1 <- alpha + n; beta.1 <- beta + sum(y)
lines(lamb,dgamma(lamb,shape = alpha.1,rate = beta.1),lty = 1,col = 'red')

alpha.1/beta.1 # Posterior Mean
pgamma(0.1,shape = alpha.1,rate = beta.1)

# ----------------- SIGNIFICANT EARTHQUAKES WORLDWIDE ------------------
# Source: The United States Geological.

# We will model the rate of earthquakes of magnitude 4.0+ in the state of California 
# during 2015. Variable to model: Waiting time (Exponential)
#   Assumptions:
#       1. Earthquake events are independents
#       2. The rate at which earthquakes occur doesn't change during the year
#       3. The earthquake hazard rate doesn't change with the time

#     Eatrhquake(1) ------------ [Time (days)] ------------> Earthquake (2)

# Prior Distribution
alpha <- 1; beta <- 30; lamb.mean <- alpha/beta; # Earthquakes per day
lamb <- 1/seq(from = 5,to = 120,by = 1);
plot(lamb,dgamma(lamb,shape = alpha,rate = beta),type = 'l',col = 'black')

# Data
dates <- as.Date(c('2015-01-4','2015-01-20','2015-01-28','2015-05-22','2015-07-21',
                   '2015-07-25','2015-08-17','2015-09-16','2015-12-30'))

# Days difference
days.diff = c()
for (j in 1:(length(dates) - 1)){
  days.diff[j] = as.numeric(dates[j + 1] - dates[j])
}
y <- days.diff
sum.y <- sum(y); n <- length(y)

# Posterior Distribution
alpha.1 <- alpha + n; beta.1 <- beta + sum.y
plot(lamb,dgamma(lamb,shape = alpha.1,rate = beta.1),type = 'l',col = 'blue')
lines(lamb,dgamma(lamb,shape = alpha,rate = beta),lty = 1,col = 'black')

lower.tail <- qgamma(0.025,shape = alpha.1,rate = beta.1)
upper.tail <- qgamma(0.975,shape = alpha.1,rate = beta.1)

# Posterior Predictive Density
y_astk <- seq(from = 0,to = 120,by = 1) # Predictive days

post_pred_exp_gamma <- function(y,alpha,beta){
  num <- alpha*beta^(alpha)
  den <- (beta + y)^(alpha + 1)
  return(num/den)
}

# Predictive Density (PPD)
PPD <- post_pred_exp_gamma(y_astk,alpha.1,beta.1)
plot(y_astk,PPD,type = 'l',col = 'blue')
