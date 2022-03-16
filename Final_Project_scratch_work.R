library(mvtnorm)
source("bsl.R")

rmodel = function(n, theta){
  #sample from the model
  return(rpois(n, theta))
}

prior_density = function(theta, log=F){
  #returns the density of the model parameter from the prior distribution
  dgamma(theta, shape=alpha, rate=beta, log=log)
}

statistic = function(X){
  #given data, this function computes a statistic
  #return(mean(X))
  #return(max(X))
  return(c(min(X), max(X)))
}

proposal = function(old_theta, new_theta=NA){
  #given the current model parameter, this function returns an update proposal
  #given the current and update, this function returns the density of the update
  bandwidth = 3
  if(is.na(new_theta)){
    return(abs(runif(1,old_theta-bandwidth,old_theta+bandwidth)))
  }
  return(1/(2*bandwidth)*ifelse(new_theta < bandwidth - old_theta, 2,1)*ifelse(abs(old_theta-new_theta)<=bandwidth,1,0))
}

lambda = 30
y = rpois(100,lambda)
sY = statistic(y)

alpha = 15
beta = 0.5
curve(dgamma(x,shape=alpha,rate=beta), 0, 60)

theta0 = 5


result = MCMC_BSL(y = y,
                  prior_density = prior_density,
                  proposal_density = function(a,b) 1,
                  rand_proposal = proposal,
                  rand_model = rmodel,
                  statistic = statistic,
                  initial_theta = theta0,
                  iterations=1000,
                  simulations=100)

view_results(result, model=list(theta=30))

plot(result$theta, type="l")
plot(matrix(unlist(result$mu), nrow=1000, ncol=2, byrow=TRUE), type="l", xlab="Min", ylab="Max", main="Mu Updates")





result_max = MCMC_BSL(Y, prior, proposal, statistic=max, n=100, N=100, iterations=1000, theta0=6.29)
result_min = MCMC_BSL(Y, prior, proposal, statistic=min, n=100, N=100, iterations=1000, theta0=6.29)
result_mean = MCMC_BSL(Y, prior, proposal, statistic=mean, n=100, N=100, iterations=1000, theta0=6.29)
plot(density(result_mean$theta[100:1000]), xlim=c(20,40), main="P(Theta|Y) for different statistics", xlab="Theta")
lines(density(result_max$theta[100:1000]), col="red")
lines(density(result_min$theta[100:1000]), col="blue")
legend(35,0.22,legend=c("Max","Mean","Min"), col=c("Red","Black","Blue"), lty=1)





#Ricker Model
rmodel_ricker = function(n, theta){
  #sample from the model
  r = exp(theta$log_r)
  sigma = theta$sigma
  phi = theta$phi

  e = rnorm(n,0,sigma)
  N = matrix(0,nrow=n)
  N[1] = 5
  Y = matrix(0,nrow=n)
  Y[1] = rpois(1,phi*N[1])
  for(i in 2:n){
    N[i] = r*N[i-1]*exp(-N[i-1]+e[i-1])
    Y[i] = rpois(1,phi*N[i])
  }
  return(Y)
}

prior_ricker = function(theta){
  #returns the density of the model parameter from the prior distribution
  return(1)
}

statistic_ricker = function(Y){
  #given data, this function computes a statistic
  stat = c(mean(Y),sum(Y==0),lm(sort(Y)~ poly(1:length(Y),2))$coefficients, arima(Y,order=c(1,0,0))$coef)
  return(as.matrix(stat,ncol=1))
}

proposal_ricker = function(old_theta, new_theta=NA){
  #given the current model parameter, this function returns an update proposal
  #given the current and update, this function returns the density of the update
  bandwidth = c(3,1,2)
  if(is.na(new_theta)){
    a = runif(1, -bandwidth[1], bandwidth[1])
    b = runif(1, -bandwidth[2], bandwidth[2])
    c = runif(1, -bandwidth[3], bandwidth[3])
    return(abs(c(old_theta+c(a,b,c))))
  }
  return(1)
}

result = MCMC_BSL(Y, rmodel_ricker, prior_ricker, proposal_ricker, statistic_ricker, n=100, N=100, iterations=100, theta=c(20,1,5))
