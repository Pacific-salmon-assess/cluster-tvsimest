compare_logbprior_func<- function(path=".", a,u){
  
  
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
             pSmax_sig=max(dat$obsSpawners)*.5
  )
  dfsip<-df
  dfsip$pSmax_sig<-max(dat$obsSpawners)*.25


  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))


  log(1/df$pSmax_mean-0.5*(df$pSmax_sig*df$pSmax_sig))
  
  f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf3<-f3$summary()

  
  f3_ip <- mod3_ip$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf3_ip<-f3_ip$summary()

  
  f3_sip <- mod3_ip$sample(data=dfsip,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf3_sip<-f3_sip$summary()


  f4 <- mod4$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4<-f4$summary()
 

  f4_ip <- mod4_ip$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4_ip<-f4_ip$summary()

  f4_sip <- mod4_ip$sample(data=dfsip,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4_sip<-f4_sip$summary()
  

  #prior values for TMB
  Smax_mean<-1/(max(df_tmb$S)*.5)
  Smax_sd<-Smax_mean
  Smax_sd_sip<-(1/(max(df_tmb$S)*.5))*.5

  ip_logb_mean<-log(Smax_mean)
  ip_logb_sd<-sqrt(log(1+ (Smax_sd)^2/ (Smax_mean)^2))
  sip_logb_sd<-sqrt(log(1+ (Smax_sd_sip)^2/ (Smax_mean)^2))

  



  ptva <- ricker_rw_TMB(data=df_tmb,tv.par='a')
  ptva_ip <- ricker_rw_TMB(data=df_tmb,tv.par='a',logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)
  ptva_sip <- ricker_rw_TMB(data=df_tmb,tv.par='a',logb_p_mean=ip_logb_mean,logb_p_sd=sip_logb_sd)

  ptvb <- ricker_rw_TMB(data=df_tmb, tv.par='b')
  ptvb_ip <- ricker_rw_TMB(data=df_tmb, tv.par='b',logb_p_mean=ip_logb_mean,logb_p_sd=ip_logb_sd)
  ptvb_sip <- ricker_rw_TMB(data=df_tmb, tv.par='b',logb_p_mean=ip_logb_mean,logb_p_sd=sip_logb_sd)
  
  #Max. prod
  
  #regime state sequence:
  
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",4),rep("MLE",4)),each=nrow(dat)),
                   model=rep(c("rwa_stan","rwa_stan_ip","rwb_stan","rwb_stan_ip","rwa_tmb","rwa_tmb_ip","rwb_tmb","rwb_tmb_ip"),each=nrow(dat)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   est=c(resf3[grep("log_a\\[",resf3$variable),"median"][[1]],
                         resf3_ip[grep("log_a\\[",resf3_ip$variable),"median"][[1]],
                         resf3_sip[grep("log_a\\[",resf3_sip$variable),"median"][[1]],
                         rep(resf4[grep("log_a",resf4$variable),"median"][[1]],nrow(dat)),
                         rep(resf4_ip[grep("log_a",resf4_ip$variable),"median"][[1]],nrow(dat)),
                         rep(resf4_sip[grep("log_a",resf4_sip$variable),"median"][[1]],nrow(dat)),
                         ptva$alpha,
                         ptva_ip$alpha,
                         ptva_sip$alpha,
                         rep(ptvb$alpha,nrow(dat)),
                         rep(ptvb_ip$alpha,nrow(dat)),
                         rep(ptvb_sip$alpha,nrow(dat))
                         ),
                   convergence= as.numeric( c(
                        abs(resf3[grep("log_a\\[",resf3$variable),"rhat"][[1]]-1)>.1,
                        abs(resf3_ip[grep("log_a\\[",resf3_ip$variable),"rhat"][[1]]-1)>.1,
                        abs(resf3_sip[grep("log_a\\[",resf3_sip$variable),"rhat"][[1]]-1)>.1,
                        rep(abs(resf4[grep("log_a",resf4$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                        rep(abs(resf4_ip[grep("log_a",resf4_ip$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                        rep(abs(resf4_sip[grep("log_a",resf4_sip$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                        rep(ptva$model$convergence + ptva$conv_problem,nrow(dat)),
                        rep(ptva_ip$model$convergence + ptva_ip$conv_problem,nrow(dat)),
                        rep(ptva_sip$model$convergence + ptva_sip$conv_problem,nrow(dat)),
                        rep(ptvb$model$convergence + ptvb$conv_problem,nrow(dat)),
                        rep(ptvb_ip$model$convergence + ptvb_ip$conv_problem,nrow(dat)),
                        rep(ptvb_sip$model$convergence + ptvb_sip$conv_problem,nrow(dat))
                        )))


  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  dfa$bias<- (dfa$est-dfa$sim)
 
  #Smax
 
  dfsmax<- data.frame(parameter="smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",4),rep("MCMC",4)),each=nrow(dat)),
                      model=rep(c("rwa_stan","rwa_stan_ip","rwb_stan","rwb_stan_ip","rwa_tmb","rwa_tmb_ip","rwb_tmb","rwb_tmb_ip"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(1/dat$beta,8),
                      est=c(rep(resf3[grep('S_max',resf3$variable),"median"][[1]],nrow(dat)),
                            rep(resf3_ip[grep('S_max',resf3_ip$variable),"median"][[1]],nrow(dat)),
                            rep(resf3_sip[grep('S_max',resf3_sip$variable),"median"][[1]],nrow(dat)),
                            resf4[grep('S_max',resf4$variable),"median"][[1]],
                            resf4_ip[grep('S_max',resf4_ip$variable),"median"][[1]],
                            resf4_sip[grep('S_max',resf4_sip$variable),"median"][[1]],
                            rep(ptva$Smax,nrow(dat)),
                            rep(ptva_ip$Smax,nrow(dat)),
                            rep(ptva_sip$Smax,nrow(dat)),
                            ptvb$Smax,
                            ptvb_ip$Smax,
                            ptvb_sip$Smax),
                      convergence=c(rep(abs(resf3[grep('S_max',resf3$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                            rep(abs(resf3_ip[grep('S_max',resf3_ip$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                            rep(abs(resf3_sip[grep('S_max',resf3_sip$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                            abs(resf4[grep('S_max',resf4$variable),"rhat"][[1]]-1)>.1,
                            abs(resf4_ip[grep('S_max',resf4_ip$variable),"rhat"][[1]]-1)>.1,
                            abs(resf4_sip[grep('S_max',resf4_sip$variable),"rhat"][[1]]-1)>.1,                   
                            rep(ptva$model$convergence + ptva$conv_problem,nrow(dat)),
                            rep(ptva_ip$model$convergence + ptva_ip$conv_problem,nrow(dat)),
                            rep(ptva_sip$model$convergence + ptva_sip$conv_problem,nrow(dat)),
                            rep(ptvb$model$convergence + ptvb$conv_problem,nrow(dat)),
                            rep(ptvb_ip$model$convergence + ptvb_ip$conv_problem,nrow(dat)),
                            rep(ptvb_sip$model$convergence + ptvb_sip$conv_problem,nrow(dat))
                            ))
            
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  dfsmax$bias<- (dfsmax$est-dfsmax$sim)
  
  
  dff<-rbind(dfa,dfsmax)
  
  return(dff)
  
}








smsy_logbprior_func<- function(path=".", a,u){
  
  
  allsimest <- list()
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
             pSmax_sig=max(dat$obsSpawners)*.5
  )


  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))


  
  #prior values for TMB
  ip_logb_mean<-log(1/(max(df_tmb$S)*.5))
  ip_logb_mean_low<-log(1/(max(df_tmb$S)*.35))
  ip_logb_sd<-sqrt(log(1+1))


  
  dff<-list(density_true_val=dnorm(log(1/unique(dat$sMSY)), ip_logb_mean, ip_logb_sd),
  density_true_val_lowpmean=dnorm(log(1/unique(dat$sMSY)), ip_logb_mean_low, ip_logb_sd),
  smsy=unique(dat$sMSY),
  maxS=max(df_tmb$S))
  
  return(dff)
  
}
