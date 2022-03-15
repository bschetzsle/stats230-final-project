library(mvtnorm)

rmodel = function(n, theta){
  #sample from the model
  return(rpois(n, theta))
}

prior = function(theta){
  #returns the density of the model parameter from the prior distribution
  dgamma(theta, shape=alpha, rate=beta)
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





MCMC_BSL = function(Y, rmodel, prior, proposal, statistic, n=100, N=100, iterations, theta){
  #this function will generate samples from the posterior P(theta|Y) without computing P(Y|theta)
  sY = statistic(Y)

  #generate a structure to keep track of the samples
  result = matrix(0, nrow=iterations, ncol=3)
  result = list(theta = numeric(), mu = list(), Sigma = list())

  #simulate x from starting theta0 and calculate statistic
  s = replicate(N, statistic(rmodel(n, theta)))
  mu = apply(s,1,sum)/n
  Sigma = (s-mu) %*% t(s-mu) / (n-1)

  for(i in 1:iterations){
    new_theta = proposal(theta)
    s = replicate(N, statistic(rmodel(n, new_theta)))
    new_mu = apply(s,1,sum)/n
    new_Sigma = (s-mu) %*% t(s-mu) / (n-1)
    r = min(1, (dmvnorm(sY, mean=new_mu, sigma=new_Sigma) * prior(new_theta) * proposal(new_theta, theta)) /
      (dmvnorm(sY, mean=mu, sigma=Sigma) * prior(theta) * proposal(theta, new_theta)) )
    if(is.na(r) || runif(1) < r){
      theta = result$theta[i] = new_theta
      mu = result$mu[[i]] = new_mu
      Sigma = result$Sigma[[i]] = new_Sigma
    }
    else{
      result$theta[i] = theta
      result$mu[[i]] = mu
      result$Sigma[[i]] = Sigma
    }
  }
  return(result)
}


lambda = 30
Y = rpois(100,lambda)
sY = statistic(Y)

alpha = 15
beta = 0.5
curve(dgamma(x,shape=alpha,rate=beta), 0, 60)

theta0 = 5


result = MCMC_BSL(Y, prior, proposal, statistic, n=100, N=100, iterations=1000, theta0=6.29)

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
  r = theta[1]
  sigma = theta[2]
  phi = theta[3]

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
