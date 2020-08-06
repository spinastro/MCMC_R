library(MASS)
library(extraDistr)
library(svMisc)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(truncnorm)
library(tmvtnorm)

set.seed(123)

################################################# likelihood
likelihood = function(param, low_lim, up_lim, obs, obs_err){
  mu1_c = param[1]
  sd1_c = param[2]
  mu2_c = param[3]
  sd2_c = param[4]
  mu3_c = param[5]
  sd3_c = param[6]
  r12_c= param[7]
  r13_c= param[8]
  r23_c= param[9]
  cov_c <-  matrix(c(sd1_c**2,r12_c,r13_c,r12_c,sd2_c**2,r23_c,r13_c,r23_c,sd3_c**2), nrow = 3, ncol = 3)
  
  mu1_f = param[10]
  sd1_f = param[11]
  mu2_f = param[12]
  sd2_f = param[13]
  mu3_f = param[14]
  sd3_f = param[15]
  r12_f= param[16]
  r13_f= param[17]
  r23_f= param[18]
  cov_f <-  matrix(c(sd1_f**2,r12_f,r13_f,r12_f,sd2_f**2,r23_f,r13_f,r23_f,sd3_f**2), nrow = 3, ncol = 3)
  
  f = param[19]
  
  y1 = rnorm(dim(obs)[1], mean = obs[,1], sd= obs_err[,1])
  y2 = rnorm(dim(obs)[1], mean = obs[,2], sd= obs_err[,2])
  y3 = rnorm(dim(obs)[1], mean = obs[,3], sd= obs_err[,3])
  Y=cbind(y1,y2,y3)
  
  if ((f > 0.001) & (f < 1)) {
    singlelikelihoods = f*(dmvnorm(Y, c(mu1_c,mu2_c,mu3_c), cov_c)/pmvnorm(lower=low_lim, upper=up_lim, mean= c(mu1_c,mu2_c,mu3_c), sigma=cov_c)) + (1.-f)*(dmvnorm(Y, c(mu1_f,mu2_f,mu3_f), cov_f)/pmvnorm(lower=low_lim, upper=up_lim, mean=c(mu1_f,mu2_f,mu3_f), sigma=cov_f))
  } else {
    singlelikelihoods = 0
  }
  #if (is.nan(singlelikelihoods[1])) {
  #  singlelikelihoods = 0  
  #}
  sumll = sum(log10(singlelikelihoods))
  return(sumll)
}

################################################# priors
prior = function(param, starting_values, prior_scales){
  p_mu1_c = param[1]
  p_sd1_c = param[2]
  p_mu2_c = param[3]
  p_sd2_c = param[4]
  p_mu3_c = param[5]
  p_sd3_c = param[6]
  p_r12_c = param[7]
  p_r13_c = param[8]
  p_r23_c = param[9]
  mu1_c_prior = dnorm(p_mu1_c, mean = starting_values[1], sd = prior_scales[1], log = T)
  sd1_c_prior = log10(dtruncnorm(p_sd1_c, a=0, b=Inf, mean = starting_values[2], sd = prior_scales[2]))
  mu2_c_prior = dnorm(p_mu2_c, mean = starting_values[3], sd = prior_scales[3], log = T)
  sd2_c_prior = log10(dtruncnorm(p_sd2_c, a=0, b=Inf, mean = starting_values[4], sd = prior_scales[4]))
  mu3_c_prior = dnorm(p_mu3_c, mean = starting_values[5], sd = prior_scales[5], log = T)
  sd3_c_prior = log10(dtruncnorm(p_sd3_c, a=0, b=Inf, mean = starting_values[6], sd = prior_scales[6]))
  r12_c_prior = dnorm(p_r12_c, mean = starting_values[7], sd = prior_scales[7], log = T)
  r13_c_prior = dnorm(p_r13_c, mean = starting_values[8], sd = prior_scales[8], log = T)
  r23_c_prior = dnorm(p_r23_c, mean = starting_values[9], sd = prior_scales[9], log = T)
  
  priors_c = c(mu1_c_prior,mu2_c_prior,mu3_c_prior,sd1_c_prior,sd2_c_prior,sd3_c_prior,r12_c_prior,r13_c_prior,r23_c_prior)
  
  p_mu1_f = param[10]
  p_sd1_f = param[11]
  p_mu2_f = param[12]
  p_sd2_f = param[13]
  p_mu3_f = param[14]
  p_sd3_f = param[15]
  p_r12_f = param[16]
  p_r13_f = param[17]
  p_r23_f = param[18]
  mu1_f_prior = dnorm(p_mu1_f, mean = starting_values[10], sd = prior_scales[10], log = T)
  sd1_f_prior = dhcauchy(p_sd1_f, sigma = prior_scales[11], log = T)
  mu2_f_prior = dnorm(p_mu2_f, mean = starting_values[12], sd = prior_scales[12], log = T)
  sd2_f_prior = dhcauchy(p_sd2_f, sigma = prior_scales[13], log = T)
  mu3_f_prior = dnorm(p_mu3_f, mean = starting_values[14], sd = prior_scales[14], log = T)
  sd3_f_prior = dhcauchy(p_sd3_f, sigma = prior_scales[15], log = T)
  r12_f_prior = dnorm(p_r12_f, mean = starting_values[16], sd = prior_scales[16], log = T)
  r13_f_prior = dnorm(p_r13_f, mean = starting_values[17], sd = prior_scales[17], log = T)
  r23_f_prior = dnorm(p_r23_f, mean = starting_values[18], sd = prior_scales[18], log = T)
  
  priors_f = c(mu1_f_prior,mu2_f_prior,mu3_f_prior,sd1_f_prior,sd2_f_prior,sd3_f_prior,r12_f_prior,r13_f_prior,r23_f_prior)
  
  p_f = param[19]
  f_prior = dnorm(p_f, mean = starting_values[19], sd = prior_scales[19], log = T)
  
  return(sum(c(priors_c,priors_f,f_prior)))
}

################################################# proposal function
proposalfunction = function(param,acc,prev_sig,step_scale){
  step = rnorm(19, mean = 0, sd= step_scale)
  if (acc == 100) {
    step[step_scale != 0] = step[step_scale != 0]*prev_sig[step_scale != 0]/sign(step[step_scale != 0])
  }
  new_par = param + step
  return(cbind(param + step, sign(step)))
}

################################################# run_metropolis_MCMC
run_metropolis_MCMC = function(iterations=1000,startvalues=rep(0,19),prior_scales=rep(1,19), step_scale=c(rep(0.05,6),rep(0.01,3),rep(0.05,6),rep(0.01,3),0.01), low_lim=c(-Inf,-Inf,-Inf), up_lim=c(Inf,Inf,Inf), burnIn=500, obs, obs_err){
  print("#######################")
  print("Start iterations")
  print("#######################")
  print("")
  
  chain = array(dim = c(iterations+1,19))
  step_values = array(dim = c(iterations+1,19))
  acceptance_probab = array(dim = iterations+1)
  acceptance_rate = array(dim = iterations+1)
  
  chain[1,] = startvalues
  step_values[1,] = step_scale
  acceptance_probab[1] = 1
  acceptance_rate[1] = 1
  
  N_acc = 0
  acc = 0
  prev_sig=c(1:19)
  for (i in 1:iterations){
    #progress(as.integer(i*100/iterations))
    #Sys.sleep(0.01)
    out_proposal = proposalfunction(chain[i,],acc,prev_sig,step_scale)
    proposal = out_proposal[1:19]
    prev_sig = out_proposal[20:38]
    step_values[i+1,] = step_scale
    acceptance_probab[i] = exp(likelihood(proposal,low_lim, up_lim, obs, obs_err) + prior(proposal,startvalues, prior_scales) - likelihood(chain[i,],low_lim, up_lim, obs, obs_err)- prior(chain[i,],startvalues, prior_scales))
    if (is.nan(acceptance_probab[i])) acceptance_probab[i] = 0
    print(i)
    print(acceptance_probab[i])
    if (runif(1) < acceptance_probab[i]) {
      chain[i+1,] = proposal
      N_acc = N_acc + 1
      #acc=1
    }else{
      chain[i+1,] = chain[i,]
      acc=0
    }
    acceptance_rate[i] = N_acc/i
    print(chain[i+1,])
    colnames(step_values) <- c("mu_1_c", "sd_1_c", "mu_2_c", "sd_2_c", "mu_3_c", "sd_3_c", "r12_c", "r13_c", "r23_c", "mu_1_f", "sd_1_f", "mu_2_f", "sd_2_f", "mu_3_f", "sd_3_f", "r12_f", "r13_f", "r23_f", "f")
    colnames(chain) <- c("mu_1_c", "sd_1_c", "mu_2_c", "sd_2_c", "mu_3_c", "sd_3_c", "r12_c", "r13_c", "r23_c", "mu_1_f", "sd_1_f", "mu_2_f", "sd_2_f", "mu_3_f", "sd_3_f", "r12_f", "r13_f", "r23_f", "f")
  }
  
  print("#######################")
  print("End iterations")
  print("#######################")
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
  print(sprintf("Acceptance: %f", acceptance))
  print("#######################")
  print("Quartiles 5% 50% 95%")
  print(sprintf("mu_1_c: %f", quantile(chain[-(1:burnIn),1], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_1_c: %f", quantile(chain[-(1:burnIn),2], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("mu_2_c: %f", quantile(chain[-(1:burnIn),3], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_2_c: %f", quantile(chain[-(1:burnIn),4], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("mu_3_c: %f", quantile(chain[-(1:burnIn),5], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_3_c: %f", quantile(chain[-(1:burnIn),6], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r12_c: %f", quantile(chain[-(1:burnIn),7], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r13_c: %f", quantile(chain[-(1:burnIn),8], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r23_c: %f", quantile(chain[-(1:burnIn),9], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("mu_1_f: %f", quantile(chain[-(1:burnIn),10], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_1_f: %f", quantile(chain[-(1:burnIn),11], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("mu_2_f: %f", quantile(chain[-(1:burnIn),12], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_2_f: %f", quantile(chain[-(1:burnIn),13], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("mu_3_f: %f", quantile(chain[-(1:burnIn),14], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("sd_3_f: %f", quantile(chain[-(1:burnIn),15], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r12_f: %f", quantile(chain[-(1:burnIn),16], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r13_f: %f", quantile(chain[-(1:burnIn),17], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("r23_f: %f", quantile(chain[-(1:burnIn),18], probs = c(0.05, 0.5, 0.95))))
  print(sprintf("f: %f", quantile(chain[-(1:burnIn),19], probs = c(0.05, 0.5, 0.95))))
  
  acceptance_probab[i+1] = acceptance_probab[i]
  acceptance_rate[i+1] = acceptance_rate[i]
  chain = cbind(chain,acceptance_probab,acceptance_rate)
  #return(chain[-(1:burnIn),])
  return(list(chain=chain,step_values=step_values))
}


posterior_sample = function(n,param, low_lim, up_lim) {
  mu1_c = param[1]
  sd1_c = param[2]
  mu2_c = param[3]
  sd2_c = param[4]
  mu3_c = param[5]
  sd3_c = param[6]
  r12_c= param[7]
  r13_c= param[8]
  r23_c= param[9]
  cov_c <-  matrix(c(sd1_c**2,r12_c,r13_c,r12_c,sd2_c**2,r23_c,r13_c,r23_c,sd3_c**2), nrow = 3, ncol = 3)
  
  mu1_f = param[10]
  sd1_f = param[11]
  mu2_f = param[12]
  sd2_f = param[13]
  mu3_f = param[14]
  sd3_f = param[15]
  r12_f= param[16]
  r13_f= param[17]
  r23_f= param[18]
  cov_f <-  matrix(c(sd1_f**2,r12_f,r13_f,r12_f,sd2_f**2,r23_f,r13_f,r23_f,sd3_f**2), nrow = 3, ncol = 3)
  
  f = param[19]
  
  y_c = rtmvnorm(as.integer(n*f), mean=c(mu1_c,mu2_c,mu3_c), sigma=cov_c, lower=low_lim, upper=up_lim)
  y_f = rtmvnorm(as.integer(n*(1-f)), mean=c(mu1_f,mu2_f,mu3_f), sigma=cov_f, lower=low_lim, upper=up_lim)
  y <- matrix(c(y_c[,1],y_f[,1],y_c[,2],y_f[,2],y_c[,3],y_f[,3]), nrow=dim(y_c)[1]+dim(y_f)[1], ncol=3, byrow=F)
  return(y)
}


