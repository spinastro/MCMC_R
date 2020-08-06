library(MASS)
library(extraDistr)

sampleSize <- 3000

true_mu <- c(10.2, 7.4)
true_sd <- c(1.5, 3.5)
r12 <- 0.2
cov_c <-  matrix(c(true_sd[1]**2,r12*true_sd[1]*true_sd[2],r12*true_sd[1]*true_sd[2],true_sd[2]**2), nrow = 2, ncol = 2)
y <- rmvnorm(n=sampleSize, true_mu, cov_c)

plot(y)

likelihood = function(param){
  mu1 = param[1]
  sd1 = param[2]
  mu2 = param[3]
  sd2 = param[4]
  r12= param[5]
  cov <-  matrix(c(sd1**2,r12*sd1*sd2,r12*sd1*sd2,sd2**2), nrow = 2, ncol = 2)
  
  singlelikelihoods = dmvnorm(y, c(mu1,mu2), cov, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

prior = function(param){
  p_mu1 = param[1]
  p_sd1 = param[2]
  p_mu2 = param[3]
  p_sd2 = param[4]
  p_r12 = param[5]
  mu1_prior = dnorm(p_mu1, mean = c(20,20), sd = c(3,3), log = T)
  sd1_prior = dhcauchy(p_sd1, sigma = c(1,1), log = T)
  mu2_prior = dnorm(p_mu2, mean = c(20,20), sd = c(3,3), log = T)
  sd2_prior = dhcauchy(p_sd2, sigma = c(1,1), log = T)
  r12_prior = dhcauchy(p_r12, sigma = 1, log = T)
  return(sum(c(mu1_prior,mu2_prior,sd1_prior,sd2_prior,r12_prior)))
}

proposalfunction = function(param){
  return(rnorm(5, mean = param, sd= c(0.1,0.1,0.1,0.1,0.1)))
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,5))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(0,2,0,2,0.1)
niter = 10000
chain = run_metropolis_MCMC(startvalue, niter)

burnIn = niter/2
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

print(tail(chain, n=1))
print(acceptance)

plot(chain[-(1:burnIn),1], type = "l")
plot(chain[,1], type = "l")

