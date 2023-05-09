emp_lfo<- function(u){
  #stocks_f
  #data_f
  dat <- data_f[data_f$stock.id2==stocks_f$stock.id2[u],]
  dat <- dat[complete.cases(dat$spawners),]

  df <- list(by=dat$broodyear,
             S=dat$spawners,
             R=dat$recruits,
             logRS=dat$klogR_S,
             L=max(dat$broodyear)-min(dat$broodyear)+1,
             ii=as.numeric(as.factor(dat$broodyear)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=c(1,1)
  )
  
  #LFO cross-validation
  #model 1 - static Ricker
  lfostatic<- samEst::stan_lfo_cv(mod=mod1lfo,type='static',df=df,L=10)
  #model 2 - static autocorrelated Ricker
  lfoac<- samEst::stan_lfo_cv(mod=mod2lfo,type='static',df=df,L=10)
  #model 3 - dynamic productivity Ricker
  lfoalpha<- samEst::stan_lfo_cv(mod=mod3lfo,type='tv',df=df,L=10)
  #model 4 - dynamic capacity Ricker
  lfobeta<- samEst::stan_lfo_cv(mod=mod4lfo,type='tv',df=df,L=10)
  #model 5 - dynamic productivity & capacity Ricker
  lfoalphabeta<- samEst::stan_lfo_cv(mod=mod5lfo,type='tv',df=df,L=10)
  #model 6 - productivity regime shift - 2 regimes
  lfohmma<- samEst::stan_lfo_cv(mod=mod6lfo,type='regime',df=df,L=10,K=2)
  #model 7 - capacity regime shift
  lfohmmb<- samEst::stan_lfo_cv(mod=mod7lfo,type='regime',df=df,L=10,K=2)
  #model 8 - productivity and capacity regime shift
  lfohmm<- samEst::stan_lfo_cv(mod=mod8lfo,type='regime',df=df,L=10,K=2)
  
  dflfo<- data.frame(parameter="LFO",
                     stock.id=u,
                     stock=stocks_f$stock.name[u] ,
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "rwa_last3","rwb_last3","rwab_last3",
                             "rwa_last5","rwb_last5","rwab_last5",
                             "hmma", "hmmb","hmmab",
                             "hmma_last3", "hmmb_last3","hmmab_last3",
                             "hmma_last5", "hmmb_last5","hmmab_last5"),
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
                     ))

  return(dflfo)
  
}