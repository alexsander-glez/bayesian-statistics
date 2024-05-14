# ----------------- Metropolis-Hastings algorithm -----------------
# For a loaded coin or fair coin, depending the number of heads (y = 2)
# in 5 flip coins

binom <- function(y,theta){
  n <- 5;
  if (theta == 1){
    f_theta = 0.6; p = 0.7;
  }
  else {f_theta = 0.4; p = 0.5}
  return(p^y*(1 - p)^(5 - y)*f_theta)
}

n <- 10000
theta <- numeric(n)
# Start with an initial value $\theta$ = 1 {loaded}
theta[1] = 1
# for i = 1,...,m
for (i in 2:n){
  # a) Propose a candidate $\theta^*$ (to be other state as $\theta_{i-1}$)
  if (theta[i - 1] == 1) {thet_astk <- 0} # $theta^*$ = 0 {fair}
  else {thet_astk <- 1}
  
  # b) Calculate alpha ratio
  g_theta_atsk <- binom(y = 2,theta = thet_astk)
  g_theta_1 <- binom(y = 2,theta = theta[i - 1])
  alph <- g_theta_atsk/g_theta_1
  
  # c) Accept proposals or reject
  if (alph > runif(1)){
    theta[i] <- thet_astk # Accept with probability alpha
  }
  else{theta[i] = theta[i - 1]}
}

mean(theta == 0) # {$\theta$ = fair}
mean(theta == 1) # {$\theta$ = loaded}
