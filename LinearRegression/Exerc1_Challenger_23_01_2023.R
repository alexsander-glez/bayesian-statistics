# ----------------- LINEAR REGRESSION ----------------
### ------------ 1986 Challenger disaster and O-rings -----

# On January 28, 1986, a routine launch was anticipated for the 
# Challenger space shuttle. Seventy-three seconds into the flight, 
# disaster happened: the shuttle broke apart, killing all seven 
# crew members on board. An investigation into the cause of the disaster 
# focused on a critical seal called an O-ring, and it is believed 
# that damage to these O-rings during a shuttle launch may be related 
# to the ambient temperature during the launch. 

# The table below summarizes observational data on O-rings for 23 shuttle 
# missions, where the mission order is based on the temperature at the time 
# of the launch.

# mission -> Shuttle mission number.
# temperature -> Temperature, in Fahrenheit.
# damaged -> Number of damaged O-rings (out of 6).
# undamaged -> Number of undamaged O-rings (out of 6). ###

origin <- read.table('https://www.openintro.org/data/tab-delimited/orings.txt',header = TRUE)
attach(origin)

# Linear regression model (lm(), with non-informative prior)
origin.lm <- lm(damaged ~ temperature)
summary(origin.lm)
length(temperature)

# Add fitted line to scatterplot
lines(temperature,fitted(origin.lm))
# 95% Posterior interval for the slope
-0.09786 - 0.02574*qt(0.975,21)
-0.09786 + 0.02574*qt(0.975,21)
confint.default(object = origin.lm)
# Note that these are the same as the frequentist confidence intervals
pairs(origin)
# The challenger launch was at 31 degrees Farenheit
# How much O-ring damage would we predict?
# y-hat
T_pred <- 31; p <- coef(origin.lm)
(y_hat <- p[1] + p[2]*T_pred)

# Posterior Prediction Interval (same as frequentist)
predict(origin.lm,data.frame(temperature = 31),interval = 'predict')
summary(origin.lm)
Ser <- 0.8521; dgf <- 21; n <- length(temperature)
y_hat - Ser*qt(0.975,dgf)*sqrt(1 + 1/n + (T_pred - mean(temperature))^2/((n-1)*var(temperature)))
y_hat + Ser*qt(0.975,dgf)*sqrt(1 + 1/n + (T_pred - mean(temperature))^2/((n-1)*var(temperature)))

# Posterior Probability that damage index is greater than zero
scle <- Ser*sqrt(1 + 1/n + (T_pred - mean(temperature))^2/((n-1)*var(temperature)))
x <- (0 - y_hat)/scle
1 - pt(x,df = 21)


# Galton's seminal data on predicting the height of children from the
# heights of the parents, all in inches
heights <- read.table('https://raw.githubusercontent.com/gupta24789/Machine-Learning-Datasets/master/Galton.txt',header = TRUE)
attach(heights)
names(heights)
pairs(heights)
plot(heights)

summary(lm(Height ~ Father + Mother + Gender + Kids))
summary(lm(Height ~ Father + Mother + Gender))
heights.lm <- lm(Height ~ Father + Mother + Gender, data = heights)

# Each extra inch taller a father is correlated with 0.4 inch extra height in the child
# Each extra inch taller a mother is correlated with 0.3 inch extra height in the child
# A male child is on average 5.3 inches taller than a female child

# 95% Posterior interval for the difference height by gender
5.226 - 0.144*qt(0.975,df = 894)
5.226 + 0.144*qt(0.975,df = 894)

# Posterior prediction interval (same as frequentist)
predict(heights.lm,data.frame(Father = 68,Mother = 64,Gender = 'M'),interval = 'predict')
predict(heights.lm,data.frame(Father = 68,Mother = 64,Gender = 'F'),interval = 'predict')
