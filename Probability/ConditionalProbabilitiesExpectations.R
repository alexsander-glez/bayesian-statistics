## ---- Conditional probabilities and expectations -------
# MIT Course: High-Dimensional Data Analysis
# url: http://genomicsclass.github.io/book/pages/conditional_expectation.html
# Sudent: González Bautista Alejandro

# Prediction problems can be divided into categorical and continuous outcomes.
## Regression in the context of prediction
# We use the son and father height example to illustrate how
# regression can be interpreted as a machine learning technique
if(!require('UsingR')){
  install.packages('UsingR')
  library(UsingR)
}
library(rafalib)
mypar(1,1)
data('father.son')
x <- round(father.son$fheight) # round to nearest inch
y <- round(father.son$sheight)
hist(y, breaks = seq(min(y), max(y)))
abline(v = mean(y), col = 'red', lwd = 2)

# We are told that the father of this randomly
# selected son has a height of 71 inches (1.25 SDs taller than
# the average). 
hist(x, breaks = seq(min(x), max(x)))
(71 - mean(x))/sd(x) # than the average
mypar(1,2)
plot(x, y, xlab = "Father's height in inches",
     ylab = "Son's height in inches",
     main = paste("correlation =", signif(cor(x,y), 2)))
abline(v= c(-0.35, 0.35) + 71, col = 'red')
hist(y[x == 71], xlab = 'Heights', nc = 8, main = "",
     xlim = range(y))
mean(y[x == 71])

# It turn out that because this data is approximated by a vivariate
# normal distribution, using calculus:
# $$f(x) = \mu_Y + \rho \frac{\sigma_Y}{\sigma_X}(X - \mu_X)$$
mypar(1,2)
plot(x, y, xlab = "Father's height in inches",
     ylab = "Son's height in inches",
     main = paste("correlation =", signif(cor(x,y), 2)))
abline(v = c(-0.35, 0.35) + 71, col = 'red')

fit <- lm(y ~ x)
abline(fit, col = 1)
hist(y[x == 71], xlab = 'Heights', nc = 8, main = "",
     xlim = range(y))
abline(v = fit$coef[1] + fit$coef[2]*71, col = 1)
