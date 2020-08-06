library(MASS)
library(extraDistr)
library(svMisc)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(reshape2)

set.seed(123)

N <- 10000
f = 0.1

true_mu_c <- c(20.2, 17.4)
true_sd_c <- c(1.5, 3.5)
true_r12_c <- 0.2
true_cov_c <-  matrix(c(true_sd_c[1]**2,true_r12_c,true_r12_c,true_sd_c[2]**2), nrow = 2, ncol = 2)

true_mu_f <- c(0.2, 1.4)
true_sd_f <- c(10.5, 13.5)
true_r12_f <- 0.05
true_cov_f <-  matrix(c(true_sd_f[1]**2,true_r12_f,true_r12_f,true_sd_f[2]**2), nrow = 2, ncol = 2)

y_c <- rmvnorm(n=as.integer(N*(f)), true_mu_c, true_cov_c)
y_f <- rmvnorm(n=as.integer(N*(1-f)), true_mu_f, true_cov_f)

y <- matrix(c(y_c[,1],y_f[,1],y_c[,2],y_f[,2]), nrow=N, ncol=2, byrow=F)

################################################# likelihood
likelihood = function(param){
  mu1_c = param[1]
  sd1_c = param[2]
  mu2_c = param[3]
  sd2_c = param[4]
  r12_c= param[5]
  cov_c <-  matrix(c(sd1_c**2,r12_c,r12_c,sd2_c**2), nrow = 2, ncol = 2)
  
  mu1_f = param[6]
  sd1_f = param[7]
  mu2_f = param[8]
  sd2_f = param[9]
  r12_f= param[10]
  cov_f <-  matrix(c(sd1_f**2,r12_f,r12_f,sd2_f**2), nrow = 2, ncol = 2)
  
  f = param[11]
  
  if ((f > 0.001) & (f < 1)) {
    singlelikelihoods = f*dmvnorm(y, c(mu1_c,mu2_c), cov_c) + (1.-f)*dmvnorm(y, c(mu1_f,mu2_f), cov_f)
  } else {
    singlelikelihoods = 0
  }
  sumll = sum(log10(singlelikelihoods))
  return(sumll)
}

################################################# priors
prior = function(param){
  p_mu1_c = param[1]
  p_sd1_c = param[2]
  p_mu2_c = param[3]
  p_sd2_c = param[4]
  p_r12_c = param[5]
  mu1_c_prior = dnorm(p_mu1_c, mean = 20, sd = 5, log = T)
  sd1_c_prior = dhcauchy(p_sd1_c, sigma = 1, log = T)
  mu2_c_prior = dnorm(p_mu2_c, mean = 20, sd = 5, log = T)
  sd2_c_prior = dhcauchy(p_sd2_c, sigma = 1, log = T)
  r12_c_prior = dnorm(p_r12_c, mean = 0, sd = 0.5, log = T)
  
  p_mu1_f = param[6]
  p_sd1_f = param[7]
  p_mu2_f = param[8]
  p_sd2_f = param[9]
  p_r12_f = param[10]
  mu1_f_prior = dnorm(p_mu1_f, mean = 0, sd = 5, log = T)
  sd1_f_prior = dhcauchy(p_sd1_f, sigma = 1, log = T)
  mu2_f_prior = dnorm(p_mu2_f, mean = 0, sd = 5, log = T)
  sd2_f_prior = dhcauchy(p_sd2_f, sigma = 1, log = T)
  r12_f_prior = dnorm(p_r12_f, mean = 0, sd = 0.5, log = T)
  
  p_f = param[11]
  f_prior = dnorm(p_f, mean = 0.5, sd = 0.1, log = T)
  
  return(sum(c(mu1_c_prior,mu2_c_prior,sd1_c_prior,sd2_c_prior,r12_c_prior,mu1_f_prior,mu2_f_prior,sd1_f_prior,sd2_f_prior,r12_f_prior,f_prior)))
}

################################################# proposal function
proposalfunction = function(param){
  step = rnorm(11, mean = 0, sd= c(0.1,0.1,0.1,0.1,0.05,0.1,0.1,0.1,0.1,0.05,0.01))
  return(param + step)
}

################################################# run_metropolis_MCMC
run_metropolis_MCMC = function(startvalue, iterations){

  chain = array(dim = c(iterations+1,11))
  chain[1,] = startvalue

  for (i in 1:iterations){
    #progress(as.integer(i*100/iterations))
    #Sys.sleep(0.01)
    proposal = proposalfunction(chain[i,])
    probab = exp(likelihood(proposal) + prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    print(i)
    print(probab)
    print(chain[i+1,])
  }
  return(chain)
}

################################################# code
startvalue = c(10,2,10,2,0.1,0,10,0,10,0.1,0.05)
niter = 10000

chain = run_metropolis_MCMC(startvalue, niter)


burnIn = niter/2
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

print(tail(chain, n=1))
print(acceptance)

plot(chain[-(1:burnIn),1], type = "l")
plot(chain[,1], type = "l")
plot(chain[,2], type = "l")
plot(chain[,3], type = "l")
plot(chain[,4], type = "l")
plot(chain[,5], type = "l")
plot(chain[,6], type = "l")
plot(chain[,7], type = "l")
plot(chain[,8], type = "l")
plot(chain[,9], type = "l")
plot(chain[,10], type = "l")
plot(chain[,11], type = "l")


print(sprintf("mu_1_c: %f", quantile(chain[-(1:burnIn),1], probs = c(0.05, 0.5, 0.95))))
print(sprintf("sd_1_c: %f", quantile(chain[-(1:burnIn),2], probs = c(0.05, 0.5, 0.95))))
print(sprintf("mu_2_c: %f", quantile(chain[-(1:burnIn),3], probs = c(0.05, 0.5, 0.95))))
print(sprintf("sd_2_c: %f", quantile(chain[-(1:burnIn),4], probs = c(0.05, 0.5, 0.95))))
print(sprintf("r12_c: %f", quantile(chain[-(1:burnIn),5], probs = c(0.05, 0.5, 0.95))))
print(sprintf("mu_1_f: %f", quantile(chain[-(1:burnIn),6], probs = c(0.05, 0.5, 0.95))))
print(sprintf("sd_1_f: %f", quantile(chain[-(1:burnIn),7], probs = c(0.05, 0.5, 0.95))))
print(sprintf("mu_2_f: %f", quantile(chain[-(1:burnIn),8], probs = c(0.05, 0.5, 0.95))))
print(sprintf("sd_2_f: %f", quantile(chain[-(1:burnIn),9], probs = c(0.05, 0.5, 0.95))))
print(sprintf("r12_f: %f", quantile(chain[-(1:burnIn),10], probs = c(0.05, 0.5, 0.95))))
print(sprintf("f: %f", quantile(chain[-(1:burnIn),11], probs = c(0.05, 0.5, 0.95))))



posterior_sample = function(n,pars) {
  mu1_c = pars[1]
  sd1_c = pars[2]
  mu2_c = pars[3]
  sd2_c = pars[4]
  r12_c= pars[5]
  cov_c <-  matrix(c(sd1_c**2,r12_c,r12_c,sd2_c**2), nrow = 2, ncol = 2)
  
  mu1_f = pars[6]
  sd1_f = pars[7]
  mu2_f = pars[8]
  sd2_f = pars[9]
  r12_f= pars[10]
  cov_f <-  matrix(c(sd1_f**2,r12_f,r12_f,sd2_f**2), nrow = 2, ncol = 2)
  
  f = pars[11]
  
  y_c = rmvnorm(as.integer(n*f), c(mu1_c,mu2_c), cov_c)
  y_f = rmvnorm(as.integer(n*(1-f)), c(mu1_f,mu2_f), cov_f)
  y <- matrix(c(y_c[,1],y_f[,1],y_c[,2],y_f[,2]), nrow=dim(y_c)[1]+dim(y_f)[1], ncol=2, byrow=TRUE)
  return(y)
}

post_pars = chain[-(1:burnIn),]

hist(y[,1], col="blue",  border="black", prob = TRUE, xlab = "x", breaks = 100)
for (i in 1:100){
  ps = posterior_sample(N, post_pars[i,])
  lines(density(ps[,1]), lwd = 1, col = "red")
}
best = c(median(chain[-(1:burnIn),1]),median(chain[-(1:burnIn),2]),median(chain[-(1:burnIn),3]),median(chain[-(1:burnIn),4]),median(chain[-(1:burnIn),5]),median(chain[-(1:burnIn),6]),median(chain[-(1:burnIn),7]),median(chain[-(1:burnIn),8]),median(chain[-(1:burnIn),9]),median(chain[-(1:burnIn),10]),median(chain[-(1:burnIn),11]))
ps = posterior_sample(N, best)
lines(density(ps[,1]), lwd = 2, col = "blue")
true = c(true_mu_c[1],true_sd_c[1],true_mu_c[2],true_sd_c[2],true_r12_c,true_mu_f[1],true_sd_f[1],true_mu_f[2],true_sd_f[2],true_r12_f,f)
ps = posterior_sample(N, true)
lines(density(ps[,1]), lwd = 2, col = "green")




