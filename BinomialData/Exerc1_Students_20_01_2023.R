# ------------------------ EXCERCISE 1 -------------------
# Suppose we are giving two students a multiple choice exam 
# with 40 questions, where each question has four choices

# We don't know how much the students have studied for this exam,
# but we think that they will do better that just guessing randomly

# 1) What are the parameters of interest?
#   - theta1: probability that the first student will answer correctly each question
#   - theta2: probability that the second student will answer correctly each question

# 2) What is our likelihood?
#   - success: theta = 1
#   - failure: theta = 0; Bernouilli (independent trials) or Binomial
#         * each question is independent
#         * probability that a student gets each question right is the same for   
#           all questions for that student (all questions have the same probability to success)

# 3) What prior should we use?
theta <- seq(from = 0.01, to = 0.99, by = 0.01)
plot(theta,dbeta(theta,1,1), type = 'l') # Uniform distribution
plot(theta, dbeta(theta,2,1), type = 'l') # effective sample size = 3
plot(theta, dbeta(theta,4,2), type = 'l') # effective sample size = 6
plot(theta, dbeta(theta,6,3), type = 'l') # effective sample size = 9
plot(theta, dbeta(theta,8,4), type = 'l') # effective sample size = 12 (more concentrated data)

# 4) What is the prior probability?
1 - pbeta(0.25,8,4) # P(theta > 0.25)
1 - pbeta(0.50,8,4) # P(theta > 0.50)
1 - pbeta(0.80,8,4) # P(theta > 0.80)

# 5) Suppose the first student gets 33 questions right. What is the posterior
# distribution for theta1?
sum.y1 <- 33
(8 + sum.y1)/(4 + 40 + 8) # Posterior mean
sum.y1/40 # MLE

# Plot the posterior first to get the right scale on the y-axis
plot(theta,dbeta(theta,8 + sum.y1,4 + 40 - sum.y1), type = 'l')
lines(theta,dbeta(theta,8,4), lty = 2)
# Plot likelihood
lines(theta,dbinom(33,40,theta), lty = 3)
# Plot scaled likelihood
lines(theta,44*dbinom(33,40,theta), lty = 3)

# Posterior probabilities
alpha.1 <- 8 + sum.y1
beta.1 <- 4 + 40 - sum.y1
1 - pbeta(0.25,alpha.1,beta.1) # P(theta1 > 0.25|Y1 = 33)
1 - pbeta(0.50,alpha.1,beta.1) # P(theta1 > 0.50|Y1 = 33)
1 - pbeta(0.80,alpha.1,beta.1) # P(theta1 > 0.80|Y1 = 33)

# Posterior credible interval
qbeta(0.025,alpha.1,beta.1)
qbeta(0.975,alpha.1,beta.1)

# 6) Suppose the second student gets 24 questions right. 
# What is the posterior distribution for theta2?
sum.y2 <- 24
alpha.2 <- 8 + sum.y2
beta.2 <- 4 + 40 - sum.y2
plot(theta,dbeta(theta,alpha.2,beta.2), type = 'l') # Posterior distribution
lines(theta,dbeta(theta,8,4), lty = 2) # Prior distribution
lines(theta,dbinom(sum.y2,40,theta),lty = 3) # Likelihood function
# Scaled likelihood
lines(theta,44*dbinom(sum.y2,40,theta),lty = 3) # Likelihood function

# Posterior mean
alpha.2/(alpha.2 + beta.2)
sum.y2/40 # MLE

# Posterior probabilities
1 - qbeta(0.25,alpha.2,beta.2) # P(theta2 > 0.25|Y2 = 24)
1 - qbeta(0.50,alpha.2,beta.2) # P(theta2 > 0.50|Y2 = 24)
1 - qbeta(0.80,alpha.2,beta.2) # P(theta2 > 0.80|Y2 = 24)

# Posterior credible interval
qbeta(0.025,alpha.2,beta.2)
qbeta(0.975,alpha.2,beta.2)

# 7) What is the posterior probability that theta1 > theta2, i.e 
# that the first student has a better chance of getting a question right 
# than the second student?
theta.1 <- rbeta(1000,alpha.1,beta.1)
theta.2 <- rbeta(1000,alpha.2,beta.2)
sum(theta.1 > theta.2)/1000
mean(theta.1 > theta.2)
