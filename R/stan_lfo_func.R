stan_lfo<- function(path=".", a,u){
  
  allsimest <- list()
  simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             logRS=log(dat$obsRecruits/dat$obsSpawners)
  )
  
  
  #options(mc.cores = 3)

  #LFO cross-validation
  #model 1 - static Ricker
  lfostatic<- stan_lfo_cv(mod=mod1lfo, type='static', df=df, L=10,
    pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  #model 2 - static autocorrelated Ricker
  lfoac<- stan_lfo_cv(mod=mod2lfo, type='static', df=df, L=10,
    pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  #model 3 - dynamic productivity Ricker
  lfoalpha<- stan_lfo_cv(mod=mod3lfo,type='tv', df=df, L=10,
    pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  #model 4 - dynamic capacity Ricker
  lfobeta<- stan_lfo_cv(mod=mod4lfo,type='tv',df=df,L=10,
    pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5,
    psig_b=max(df$S)*.5)
  #model 5 - dynamic productivity & capacity Ricker
  lfoalphabeta<- stan_lfo_cv(mod=mod5lfo,type='tv',df=df,L=10,
    pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5,
    psig_b=max(df$S)*.5)
  #model 6 - productivity regime shift - 2 regimes
  lfohmma<- stan_lfo_cv(mod=mod6lfo,type='regime',df=df,L=10,K=2,
    dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2), pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  #model 7 - capacity regime shift
  lfohmmb<- stan_lfo_cv(mod=mod7lfo,type='regime',df=df,L=10,K=2,
    dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2), pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  #model 8 - productivity and capacity regime shift
  lfohmm<- stan_lfo_cv(mod=mod8lfo,type='regime',df=df,L=10,K=2,
    dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2), pSmax_mean=max(df$S)*.5,
    pSmax_sig=max(df$S)*.5)
  
  dflfo<- data.frame(parameter="LFO",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",20),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "rwa_last3","rwb_last3","rwab_last3",
                             "rwa_last5","rwb_last5","rwab_last5",
                             "hmma", "hmmb","hmmab",
                             "hmma_last3", "hmmb_last3","hmmab_last3",
                             "hmma_last5", "hmmb_last5","hmmab_last5"),
                     by=rep(NA,20),
                     sim=rep(NA,20),
                     est=c(sum(lfostatic), 
                           sum(lfoac), 
                           sum(lfoalpha[1,]), 
                           sum(lfoalpha[2,]), 
                           sum(lfoalpha[3,]), 
                           sum(lfobeta[1,]), 
                           sum(lfobeta[2,]), 
                           sum(lfobeta[3,]), 
                           sum(lfoalphabeta[1,]), 
                           sum(lfoalphabeta[2,]), 
                           sum(lfoalphabeta[3,]),    
                           sum(lfohmma[1,]),
                           sum(lfohmma[2,]),
                           sum(lfohmma[3,]),
                           sum(lfohmmb[1,]),
                           sum(lfohmmb[2,]),
                           sum(lfohmmb[3,]),
                           sum(lfohmm[1,]),
                           sum(lfohmm[2,]),
                           sum(lfohmm[3,])
                     ),
                     convergence=rep(NA,20),
                     pbias=rep(NA,20))
  
  return(dflfo)
  
}

