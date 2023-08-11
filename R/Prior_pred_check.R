

library(ggplot2)
library(cmdstanr)
library(rstan)
library(samEst)
library(posterior)
#data

a=6
u=15
path="."
  simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=max(dat$obsSpawners)*.5,
             pSmax_sig=max(dat$obsSpawners)*2,
             psig_b=.4
  )

logbeta_pr_sig=sqrt(log(1+((1/df$pSmax_sig)*(1/df$pSmax_sig))/((1/df$pSmax_mean)*(1/df$pSmax_mean)))); #this converts sigma on the untransformed scale to a log scale
logbeta_pr=log(1/df$pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; #//convert smax prior to per capita slope - transform to log scale with bias correction

#rwb


#compile 
file4ppc=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ppc.stan")
mod4ppc=cmdstanr::cmdstan_model(file4ppc)

ppc4<- mod4ppc$sample(data=df,
                    fixed_param=TRUE)
f4_ppc<-ppc4$summary()
  
prch<-ppc4$draws("sigma")
head(prch)

bayesplot::mcmc_areas(
  ppc4$draws(colnames(draws)[grep("S_max\\[",colnames(draws))[1:10]]), 
  prob = .05,
   prob_outer = 0.9,
   point_est = "mean"
)+
coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  


bayesplot::mcmc_areas(
  ppc4$draws(colnames(draws)[grep("S_max\\[",colnames(draws))[31:40]]), 
  prob = 2/3,
   prob_outer = 0.9,
   point_est = "mean"
)+
coord_cartesian(xlim = c(70000,1000000))+ 
theme_bw(12)
  

bayesplot::mcmc_areas(
  ppc4$draws(colnames(draws)[grep("S_max\\[",colnames(draws))]), 
  prob = .05,
   prob_outer = 0.9,
   point_est = "median"
)+
#coord_cartesian(xlim = c(70000,1000000))+ 
theme_bw(12)
  

bayesplot::mcmc_areas(
  ppc4$draws(c("sigma_b","sigma")), 
  prob = 2/3,
   prob_outer = 0.9,
   point_est = "mean"
)+
#coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  



bayesplot::mcmc_areas(
  ppc4$draws("log_a"), 
  prob = .05,
   prob_outer = 0.9,
   point_est = "mean"
)+
#coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  




#rw_smax

df_smax <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=max(dat$obsSpawners)*.5,
             pSmax_sig=max(dat$obsSpawners)*.5,
             psig_b=max(dat$obsSpawners)*.25
  )
#compile 
file4smax_ppc=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_smax_ppc.stan")
mod4smax_ppc=cmdstanr::cmdstan_model(file4smax_ppc)

ppc4_smax<- mod4smax_ppc$sample(data=df_smax,
                    fixed_param=TRUE)
f4_ppc_smax<-ppc4_smax$summary()
drawsmax<-ppc4_smax$draws()  

names(drawsmax)

prch<-ppc4_smax$draws("sigma")
head(prch)

bayesplot::mcmc_areas(
  ppc4_smax$draws(),
  regex_pars = "Smax\\[", 
  prob = .05,
   prob_outer = 0.9,
   point_est = "mean"
)+
#coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  


bayesplot::mcmc_areas(
  ppc4$draws(colnames(draws)[grep("S_max\\[",colnames(draws))[31:40]]), 
  prob = 2/3,
   prob_outer = 0.9,
   point_est = "mean"
)+
coord_cartesian(xlim = c(70000,1000000))+ 
theme_bw(12)
  

bayesplot::mcmc_areas(
  ppc4$draws(colnames(draws)[grep("S_max\\[",colnames(draws))]), 
  prob = .05,
   prob_outer = 0.9,
   point_est = "median"
)+
#coord_cartesian(xlim = c(70000,1000000))+ 
theme_bw(12)
  

bayesplot::mcmc_areas(
  ppc4_smax$draws(c("sigma_b","sigma")), 
  prob = 2/3,
   prob_outer = 0.9,
   point_est = "mean"
)+
#coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  



bayesplot::mcmc_areas(
  ppc4$draws("log_a"), 
  prob = .05,
   prob_outer = 0.9,
   point_est = "mean"
)+
#coord_cartesian(xlim = c(70000,500000))+ 
theme_bw(12)
  


