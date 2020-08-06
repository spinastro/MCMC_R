library(extraDistr)

true_mu <- 10.5
true_sigma <- 3.2
sampleSize <- 3000

y <-  rnorm(n=sampleSize, mean=true_mu, sd=true_sigma)

hist(y, breaks=30)

likelihood = function(param){
  mu = param[1]
  sd = param[2]
  
  singlelikelihoods = dnorm(y, mean = mu, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

prior = function(param){
  p_mu = param[1]
  p_sd = param[2]
  mu_prior = dnorm(p_mu, mean = 20, sd = 3, log = T)
  sd_prior = dhcauchy(p_sd, sigma = 1, log = T)
  return(mu_prior + sd_prior)
}

proposalfunction = function(param){
  return(rnorm(2, mean = param, sd= c(0.1,0.1)))
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
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

startvalue = c(-4,20)
chain = run_metropolis_MCMC(startvalue, 100000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

tail(chain, n=1)
print(acceptance)

plot(chain)
