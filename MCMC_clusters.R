source('MCMC_3D_mvmg.R')
source('MCMC_3D_mvmg_robust.R')
library('FITSio')
library(astroFns)

SkyDist = function(ra1,dec1,ra2,dec2) {
  # Input coordinates in deg. Output in deg.
  ra1 <- ra1*pi/180
  dec1 <- dec1*pi/180
  ra2 <- ra2*pi/180
  dec2 <- dec2*pi/180
  a <- acos(sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2))
  a*180./pi
}

Parallax_settings = function(Gaia,Plx,s_Plx) {
  explore_plx = Gaia$parallax[((Gaia$parallax > Plx+2*s_Plx) | (Gaia$parallax < Plx-2*s_Plx)) & (Gaia$parallax <2)]
  hist_plx = hist(explore_plx,50)
  peak_plx = hist_plx$mids[hist_plx$density == max(hist_plx$density)][1]
  abline(v = peak_plx, col="red")
  
  sigma_close = 5
  sigma_intermediate = 7.
  sigma_far = 3.
  
  if (Plx > 2.) { # Nearby cluster
    print("Nearby cluster")
    L_plx = Plx - 10 * s_Plx
    U_plx = Plx + 5 * s_Plx
    mu_plx_field = peak_plx
    sd_plx_field = 2
  } else if ((Plx < 2.) & (Plx > 1)) {  # Intermediate distance cluster
    print("Intermediate distance cluster")
    L_plx = Plx - 7 * s_Plx
    U_plx = Plx + 7 * s_Plx
    mu_plx_field = peak_plx
    sd_plx_field = 1
  } else if ((Plx < 1.) & (Plx > peak_plx)) {  # Distant cluster
    print("Distant cluster")
    L_plx = Plx - 7 * s_Plx
    U_plx = Plx + 7 * s_Plx
    mu_plx_field = peak_plx
    sd_plx_field = 1
  } else if (Plx < peak_plx) { # Very distant cluster
    print("Very distant cluster")
    L_plx = Plx - 3 * s_Plx
    U_plx = Plx + 3 * s_Plx
    if (L_plx<0.2) L_plx = 0.2
    mu_plx_field = peak_plx
    sd_plx_field = peak_plx/2
  } else if (Plx < 0.3) { # Very distant cluster
    print("Very distant cluster")
    L_plx = 0.2
    U_plx = Plx + 2 * s_Plx
    mu_plx_field = peak_plx
    sd_plx_field = peak_plx/2
  }
  return(c(L_plx,U_plx,mu_plx_field,sd_plx_field,ss_mu1_f,ss_sd1_f,peak_plx))
}


pm_Settings = function(Gaia,L_plx,U_plx,pmRA,s_pmRA,pmDE,s_pmDE) {

  y_pmRA = Gaia$pmra[(Gaia$parallax > L_plx) & (Gaia$parallax < U_plx)]
  y_pmDEC = Gaia$pmdec[(Gaia$parallax > L_plx) & (Gaia$parallax < U_plx)]
  
  explore_pmRA = y_pmRA[((y_pmRA > pmRA+3*s_pmRA) | (y_pmRA < pmRA-3*s_pmRA))]
  hist_pmRA = hist(explore_pmRA,100)
  peak_pmRA = hist_pmRA$mids[hist_pmRA$density == max(hist_pmRA$density)]
  sd_pmRA = sd(explore_pmRA)
  explore_pmDEC = y_pmDEC[((y_pmDEC > pmDE+3*s_pmDE) | (y_pmDEC < pmDE-3*s_pmDE))]
  hist_pmDEC = hist(explore_pmDEC,100)
  peak_pmDEC = hist_pmDEC$mids[hist_pmDEC$density == max(hist_pmDEC$density)]
  sd_pmDEC = sd(explore_pmDEC)
  
  par(mfrow=c(3,1))
  hist(Gaia$parallax,1000, xlim=c(0.0,7))
  abline(v = Plx, col="red")
  abline(v = peak_plx, col="blue")
  hist(Gaia$pmra,5000, xlim=c(-30,50))
  abline(v = pmRA, col="red")
  abline(v = peak_pmRA, col="blue")
  abline(v = peak_pmRA + sd_pmRA, col="blue")
  abline(v = peak_pmRA - sd_pmRA, col="blue")
  hist(Gaia$pmdec,5000, xlim=c(-30,30))
  abline(v = pmDE, col="red")
  abline(v = peak_pmDEC, col="blue")
  abline(v = peak_pmDEC + sd_pmDEC, col="blue")
  abline(v = peak_pmDEC - sd_pmDEC, col="blue")
  par(mfrow=c(1,1))
  
  L_pmRA=peak_pmRA-2*sd_pmRA
  U_pmRA=peak_pmRA+2*sd_pmRA
  L_pmDEC=peak_pmDEC-2*sd_pmDEC
  U_pmDEC=peak_pmDEC+2*sd_pmDEC
  return(c(L_pmRA,U_pmRA,L_pmDEC,U_pmDEC,peak_pmRA,sd_pmRA,peak_pmDEC,sd_pmDEC))
}


GaiaDR2 <- readFITS("~/Projects/GALAH/find_members/tables/GaiaDR2_clusters.fit")


catalogs = list.files(path="~/Projects/GALAH/find_members/tables/", pattern = "_2s.csv", recursive = TRUE)
log_iterations = array(dim = c(length(catalogs),80))
cluster_params = c("Cluster","RA_ICRS","DE_ICRS","r50","pmRA","s_pmRA","e_pmRA","pmDE","s_pmDE","e_pmDE","Plx","s_Plx","e_Plx","sigma_radius")
limits = c("ll_plx","ll_pmRA","ll_pmDE","ul_plx","ul_pmRA","ul_pmDE")
startvalues_names = c("sv_mu_1_c", "sv_sd_1_c", "sv_mu_2_c", "sv_sd_2_c", "sv_mu_3_c", "sv_sd_3_c", "sv_r12_c", "sv_r13_c", "sv_r23_c", "sv_mu_1_f", "sv_sd_1_f", "sv_mu_2_f", "sv_sd_2_f", "sv_mu_3_f", "sv_sd_3_f", "sv_r12_f", "sv_r13_f", "sv_r23_f", "sv_f")
prior_scales_names = c("ps_mu_1_c", "ps_sd_1_c", "ps_mu_2_c", "ps_sd_2_c", "ps_mu_3_c", "ps_sd_3_c", "ps_r12_c", "ps_r13_c", "ps_r23_c", "ps_mu_1_f", "ps_sd_1_f", "ps_mu_2_f", "ps_sd_2_f", "ps_mu_3_f", "ps_sd_3_f", "ps_r12_f", "ps_r13_f", "ps_r23_f", "ps_f")
step_scale_names = c("ss_mu_1_c", "ss_sd_1_c", "ss_mu_2_c", "ss_sd_2_c", "ss_mu_3_c", "ss_sd_3_c", "ss_r12_c", "ss_r13_c", "ss_r23_c", "ss_mu_1_f", "ss_sd_1_f", "ss_mu_2_f", "ss_sd_2_f", "ss_mu_3_f", "ss_sd_3_f", "ss_r12_f", "ss_r13_f", "ss_r23_f", "ss_f")
colnames(log_iterations) = c(cluster_params,limits,startvalues_names,prior_scales_names,step_scale_names,"Iterations", "BurnIn", "acceptance")


i=7

for (i in 1:length(catalogs)) {
  
  ####### Find cluster's parameters
  name = unlist(strsplit(catalogs[i],"_2s"))[1]
  print(name)
  
  f = substr(GaiaDR2$DF$Cluster, 1, nchar(name)) == name
  cluster = GaiaDR2$DF$Cluster[f]
  RA_ICRS = GaiaDR2$DF$RA_ICRS[f]*1000
  DE_ICRS = GaiaDR2$DF$DE_ICRS[f]*1000
  r50 = GaiaDR2$DF$r50[f] * 1000
  pmRA = GaiaDR2$DF$pmRA[f] * 1000
  s_pmRA = GaiaDR2$DF$s_pmRA[f] * 1000
  e_pmRA = GaiaDR2$DF$e_pmRA[f] * 1000
  pmDE = GaiaDR2$DF$pmDE[f] * 1000
  s_pmDE = GaiaDR2$DF$s_pmDE[f] * 1000
  e_pmDE = GaiaDR2$DF$e_pmDE[f] * 1000
  Plx = GaiaDR2$DF$Plx[f] * 1000
  s_Plx = GaiaDR2$DF$s_Plx[f] * 1000
  e_Plx = GaiaDR2$DF$e_Plx[f] * 1000
  
  
  ####### Open Gaia database
  mydata <- read.table(paste0("~/Projects/GALAH/find_members/tables/",catalogs[i]), header=TRUE,sep=",")
  distances = SkyDist(RA_ICRS, DE_ICRS, mydata$ra, mydata$dec)
  
  lim_plx = 0.
  lim_rel_err_plx = 0.2
  lim_G = 18
  sigma_radius=2
  lim_distance = sigma_radius*r50/0.6745
  Gaia = mydata[(mydata$parallax > lim_plx) & (mydata$parallax_error/mydata$parallax < lim_rel_err_plx) & (mydata$phot_g_mean_mag < lim_G) & (distances < lim_distance),]
  
  if (dim(Gaia)[1] < 500) {
    sigma_radius=3
    lim_distance = sigma_radius*r50/0.6745
    Gaia = mydata[(mydata$parallax > lim_plx) & (mydata$parallax_error/mydata$parallax < lim_rel_err_plx) & (mydata$phot_g_mean_mag < lim_G) & (distances < lim_distance),]
  }
  
  par(mfrow=c(3,1))
  hist(Gaia$parallax,3000, xlim=c(0.0,10))
  abline(v = Plx, col="red")
  hist(Gaia$pmra,3000, xlim=c(-30,50))
  abline(v = pmRA, col="red")
  hist(Gaia$pmdec,3000, xlim=c(-30,30))
  abline(v = pmDE, col="red")
  par(mfrow=c(1,1))
  
  ####### Set upper/lower limit in Plx
  array = Parallax_settings(Gaia,Plx,s_Plx)
  L_plx= array[1]
  U_plx= array[2]
  mu_plx_field= array[3]
  sd_plx_field= array[4]
  ss_mu1_f= array[5]
  ss_sd1_f= array[6]
  peak_plx= array[7]
  
  ####### Set upper/lower limit in pmRA/pmDEC
  array = pm_Settings(Gaia,L_plx,U_plx,pmRA,s_pmRA,pmDE,s_pmDE)
  L_pmRA = array[1]
  U_pmRA = array[2]
  L_pmDEC = array[3]
  U_pmDEC = array[4]
  peak_pmRA = array[5]
  sd_pmRA = array[6]
  peak_pmDEC = array[7]
  sd_pmDEC = array[8]
  
  ####### Set lower and upper limits
  low_lim = c(L_plx,L_pmRA,L_pmDEC) # Lower limits of the truncation
  up_lim = c(U_plx,U_pmRA,U_pmDEC) # Upper limits of the truncation
  
  ####### Select 
  y_precut = cbind(Gaia$parallax,Gaia$pmra,Gaia$pmdec)
  y_precut_error = cbind(Gaia$parallax_error,Gaia$pmra_error,Gaia$pmdec_error)
  filter = (y_precut[,1] > low_lim[1]) & (y_precut[,1] < up_lim[1]) & (y_precut[,2] > low_lim[2]) & (y_precut[,2] < up_lim[2]) & (y_precut[,3] > low_lim[3]) & (y_precut[,3] < up_lim[3])
  y = y_precut[filter,,drop=FALSE]
  y_err = y_precut_error[filter,,drop=FALSE]
  y_err = cbind(rep(0,dim(y)[1]),rep(0,dim(y)[1]),rep(0,dim(y)[1]))
  
  pdf(file = paste0("../plots/alldata_",name,".pdf"))
  par(mfrow=c(3,1))
  hist(y_precut[,1],1000, xlim=c(0,Plx + 10*s_Plx))
  abline(v = Plx, col="red")
  abline(v = Plx - s_Plx, col="red")
  abline(v = Plx + s_Plx, col="red")
  abline(v = peak_plx, col="blue")
  hist(y_precut[,2],500, xlim=c(peak_pmRA - 4*sd_pmRA,peak_pmRA + 4*sd_pmRA))
  abline(v = pmRA, col="red")
  abline(v = pmRA - s_pmRA, col="red")
  abline(v = pmRA + s_pmRA, col="red")
  abline(v = peak_pmRA, col="blue")
  abline(v = peak_pmRA + sd_pmRA, col="blue")
  abline(v = peak_pmRA - sd_pmRA, col="blue")
  hist(y_precut[,3],500, xlim=c(peak_pmDEC - 4*sd_pmDEC,peak_pmDEC + 4*sd_pmDEC))
  abline(v = pmDE, col="red")
  abline(v = pmDE - s_pmDE, col="red")
  abline(v = pmDE + s_pmDE, col="red")
  abline(v = peak_pmDEC, col="blue")
  abline(v = peak_pmDEC + sd_pmDEC, col="blue")
  abline(v = peak_pmDEC - sd_pmDEC, col="blue")
  par(mfrow=c(1,1))
  dev.off()
  
  pdf(file = paste0("../plots/obs_",name,".pdf"))
  par(mfrow=c(3,1))
  hist(y[,1],100)
  abline(v = Plx, col="red")
  abline(v = peak_plx, col="blue")
  hist(y[,2],500)
  abline(v = pmRA, col="red")
  abline(v = peak_pmRA, col="blue")
  abline(v = peak_pmRA + sd_pmRA, col="blue")
  abline(v = peak_pmRA - sd_pmRA, col="blue")
  hist(y[,3],500)
  abline(v = pmDE, col="red")
  abline(v = peak_pmDEC, col="blue")
  abline(v = peak_pmDEC + sd_pmDEC, col="blue")
  abline(v = peak_pmDEC - sd_pmDEC, col="blue")
  par(mfrow=c(1,1))
  dev.off()
  
  ####### Set start values
  sv_c = c(Plx, s_Plx , pmRA , s_pmRA, pmDE, s_pmDE, 0,0,0)
  sv_f = c(mu_plx_field, sd_plx_field , peak_pmRA , sd_pmRA, peak_pmDEC, sd_pmDEC, 0,0,0)
  sv_frac = 0.3
  startvalues = c(sv_c, sv_f, sv_frac) # starting values. This also sets the mu parameters in the priors.
  
  ####### Set prior scales
  ps_c = c(0.01,0.001,0.01,0.001,0.01,0.001,0.01,0.01,0.01)
  ps_f = c(0.1,10,0.1,10,0.1,10,0.1,0.1,0.1)
  ps_frac = 0.2
  prior_scales = c(ps_c,ps_f,ps_frac) # The scale parameters in the priors.
  
  ####### Set step scales
  ss_c = c(0.01,0.001,0.01,0.001,0.01,0.001,0.01,0.01,0.01)
  ss_f = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
  ss_frac = 0.01
  step_scale = c(ss_c,ss_f,ss_frac) # Scale of each step in the parameter space. The parameter won't change if its scale is 0.
  
  ####### Start MCMC
  n_ite=100000 # Number of iterations
  burnIn = n_ite/2 # Burn-in
  #trace = run_metropolis_MCMC(iterations=n_ite, startvalues=startvalues, prior_scales=prior_scales, step_scale=step_scale, low_lim=low_lim, up_lim=up_lim, burnIn, obs=y, obs_err = y_err)
  trace = run_metropolis_MCMC_robust(iterations=n_ite, startvalues=startvalues, prior_scales=prior_scales, step_scale=step_scale, low_lim=low_lim, up_lim=up_lim, S=diag(1,19), burnIn, obs=y, obs_err = y_err, adapt=TRUE)
  acceptance = 1-mean(duplicated(trace))
    
  write.csv(trace$chain, file = paste0("~/Projects/GALAH/find_members/tables/trace_",name,".csv"), quote=F, row.names = F)
  write.csv(trace$step_values, file = paste0("~/Projects/GALAH/find_members/tables/steps_",name,".csv"), quote=F, row.names = F)
  
  log_iterations[i,] = c(name,RA_ICRS,DE_ICRS,r50,pmRA,s_pmRA,e_pmRA,pmDE,s_pmDE,e_pmDE,Plx,s_Plx,e_Plx,sigma_radius,low_lim,up_lim,startvalues,prior_scales,step_scale,n_ite,burnIn,acceptance)

}

write.csv(log_iterations, file = "~/Projects/GALAH/find_members/tables/log_iterations.csv", quote=F, row.names = F)

####################################################################


plot(trace$chain[,11], type='l')




