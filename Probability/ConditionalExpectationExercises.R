## ---------- Conditional Expectation Exercises --------------
# MIT Course: High-Dimensional Data Analysis
# url: https://learning.edx.org/course/course-v1:HarvardX+PH525.4x+2T2022/block-v1:HarvardX+PH525.4x+2T2022+type@sequential+block@20be257ec0c341309492de6a7f2e695a/block-v1:HarvardX+PH525.4x+2T2022+type@vertical+block@4c469a10b04f434b8793397006bba4fc
# Conditional Expectations
# Student: González Bautista Alejandro

# Throughout this assessment it will be useful to remember
# that when our data are 0s and 1s, probabilities and expectations
# are the same thing.
n <- 1000
y <- rbinom(n, size = 1, prob = 0.25)
# Proportion of ones Pr(Y)
sum(y == 1)/length(y)
# Expectation of Y
mean(y)

# |->| Conditional Expectation Exercises #1
n <- 10000
set.seed(1)
men <- rnorm(n, 176, 7) # Heights in centimeters
women <- rnorm(n, 162, 7) # Heights in centimeters
y <- c(rep(0,n), rep(1,n))
x <- round(c(men, women))
# Mix it up
ind <- sample(seq(along = y))
y = y[ind] # labels: men(0), women(1)
x = x[ind] # Heights data

# Treating the data generated above as the population, if
# we know someone is 176 cm tall, what it the probability
# that this person is a woman:
mean(y[x == 176])


# |->| Conditional Expectation Exercises #2
# Now make a plot of E(Y|X = x) for `x = seq(160,178)`
# using the data generated in Conditional Expectation
# Exercises #1
#
# Suppose for each height x you predict 1 (female) if
# Pr(Y|X = x) > 0.5 and 0 (male) otherwise.
# What is the larger height for you predict female?
xs <- seq(160, 178)
pr <- sapply(xs, function(x0) mean(y[x == x0]))
plot(xs,pr)
abline(h = 0.5)
abline(v = 168)
plot(xs, 1 - pr)
