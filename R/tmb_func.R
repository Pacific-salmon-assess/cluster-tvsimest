
tmb_func <- function(path=".",a, u) {
  
  
  simData <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  #compiled Bayesian models try moving this out of function
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))
  
  dirpr<-matrix(c(2,1,1,2),2,2)

  p <- ricker_TMB(data=df)
  pac <- ricker_TMB(data=df, AC=TRUE)
  ptva <- ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- ricker_rw_TMB(data=df, tv.par='both')
  phmma <- ricker_hmm_TMB(data=df, tv.par='a', dirichlet_prior=dirpr)
  phmmb <- ricker_hmm_TMB(data=df, tv.par='b', dirichlet_prior=dirpr)
  phmm  <- ricker_hmm_TMB(data=df, tv.par='both', dirichlet_prior=dirpr)

  

  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",8)),each=nrow(df)),
              model=rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=nrow(df)),
              by=rep(dat$year,8),
              sim=rep(dat$alpha,8),
              median=NA,
              mode=c(rep(p$alpha,nrow(df)),
                    rep(pac$alpha,nrow(df)),
                    ptva$alpha,
                    rep(ptvb$alpha,nrow(df)),
                    ptvab$alpha,
                    phmma$alpha[phmma$regime],
                    rep(phmmb$alpha,nrow(df)),
                    phmm$alpha[phmm$regime]), 
              convergence=c(rep(c(p$model$convergence + p$conv_problem,
                    pac$model$convergence + pac$conv_problem,
                    ptva$model$convergence + ptva$conv_problem,
                    ptvb$model$convergence + ptvb$conv_problem,
                    ptvab$model$convergence + ptvab$conv_problem,
                    phmma$model$convergence + phmma$conv_problem,
                    phmmb$model$convergence + phmmb$conv_problem,
                    phmm$model$convergence + phmm$conv_problem
                    ),each=nrow(df))))
                    
  dfa$pbias<- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias<- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",8)),each=nrow(df)),
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma","hmmb","hmmab"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(1/dat$beta,8),
      median=NA,
      mode=c(rep(p$Smax,nrow(df)),
        rep(pac$Smax,nrow(df)),
        rep(ptva$Smax,nrow(df)),
        ptvb$Smax,
        ptvab$Smax,
        rep(phmma$Smax,nrow(df)),
        phmmb$Smax[phmmb$regime],
        phmm$Smax[phmm$regime]),
      convergence=rep(c(p$model$convergence + p$conv_problem,
        pac$model$convergence + pac$conv_problem,
        ptva$model$convergence + ptva$conv_problem,
        ptvb$model$convergence + ptvb$conv_problem,
        ptvab$model$convergence + ptvab$conv_problem,
        phmma$model$convergence + phmma$conv_problem,
        phmmb$model$convergence + phmmb$conv_problem,
        phmm$model$convergence + phmm$conv_problem),each=nrow(df)))
      
    dfsmax$pbias<- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias<- (dfsmax$mode-dfsmax$sim)
       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma","hmmb","hmmab"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(dat$sigma,8),
      median=NA,
      mode=c(rep(p$sig,nrow(df)),
        rep(pac$sig,nrow(df)),
        rep(ptva$sig,nrow(df)),
        rep(ptvb$sig,nrow(df)),
        rep(ptvab$sig,nrow(df)),
        rep(phmma$sigma,nrow(df)),
        rep(phmmb$sigma,nrow(df)),
        rep(phmm$sigma,nrow(df))),
      
      convergence=rep(c(p$model$convergence + p$conv_problem,
        pac$model$convergence + pac$conv_problem,
        ptva$model$convergence + ptva$conv_problem,
        ptvb$model$convergence + ptvb$conv_problem,
        ptvab$model$convergence + ptvab$conv_problem,
        phmma$model$convergence + phmma$conv_problem,
        phmmb$model$convergence + phmmb$conv_problem,
        phmm$model$convergence + phmm$conv_problem),each=nrow(df)))
    
    dfsig$pbias<- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias<- (dfsig$mode-dfsig$sim)

    #sigma_a
    dfsiga<- data.frame(parameter="sigma_a",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=c("rwa","rwab"),
      by=NA,
      sim=NA,
      median=NA,
      mode=c(ptva$siga,ptvab$siga),
      convergence=c( ptva$model$convergence + ptva$conv_problem,
        
        ptvab$model$convergence + ptvab$conv_problem))
    
    dfsiga$pbias<- ((dfsiga$mode-dfsiga$sim)/dfsiga$sim)*100
    dfsiga$bias<- (dfsiga$mode-dfsiga$sim)

    #sigma_b
    dfsigb<- data.frame(parameter="sigma_b",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=c("rwb","rwab"),
      by=NA,
      sim=NA,
      median=NA,
      mode=c(ptvb$sigb,ptvab$sigb),
      convergence=c( ptvb$model$convergence + ptvb$conv_problem,        
        ptvab$model$convergence + ptvab$conv_problem))
    
    dfsigb$pbias<- ((dfsigb$mode-dfsigb$sim)/dfsigb$sim)*100
    dfsigb$bias<- (dfsigb$mode-dfsigb$sim)
              
    #Smsy
    smsysim<-smsyCalc(dat$alpha,dat$beta)
  
    dfsmsy<- data.frame(parameter="smsy",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",8)),each=nrow(df)),
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma","hmmb","hmmab"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(smsysim,8),
      median=NA,
      mode=c(rep(p$Smsy,nrow(df)),
            rep(pac$Smsy,nrow(df)),
            ptva$Smsy,
            ptvb$Smsy,
            ptvab$Smsy,
            phmma$Smsy[phmma$regime],
            phmmb$Smsy[phmmb$regime],
            phmm$Smsy[phmm$regime]),    
      convergence=rep(c(p$model$convergence + p$conv_problem,
        pac$model$convergence + pac$conv_problem,
        ptva$model$convergence + ptva$conv_problem,
        ptvb$model$convergence + ptvb$conv_problem,
        ptvab$model$convergence + ptvab$conv_problem,
        phmma$model$convergence + phmma$conv_problem,
        phmmb$model$convergence + phmmb$conv_problem,
        phmm$model$convergence + phmm$conv_problem),each=nrow(df))) 
  
  dfsmsy$pbias<- ((dfsmsy$mode-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (dfsmsy$mode-dfsmsy$sim)

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",8)),each=nrow(df)),
    model=rep(c("simple",
      "autocorr",
      "rwa","rwb","rwab",
      "hmma","hmmb","hmmab"),each=nrow(df)),
    by=rep(dat$year,8),
    sim=rep(unlist(mapply(sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),8),
    median=NA,
    mode=c(unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="simple"&dfa$method=="MLE"],
            Smsy=dfsmsy$mode[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
            b=1/dfsmax$mode[dfsmax$model=="simple"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="autocorr"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="autocorr"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwa"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"]))),
     convergence=rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptvab$conv_problem,
      ptvb$model$convergence + ptvab$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem),
      each=nrow(df)))
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
  #1-gsl::lambert_W0(exp(1 - ptva$alpha))) /ptva$beta

    dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",8)),each=nrow(df)),
    model=rep(c("simple",
      "autocorr",
      "rwa","rwb","rwab",
      "hmma","hmmb","hmmab"),each=nrow(df)),
    by=rep(dat$year,8),
    sim=rep(umsyCalc(dat$alpha),8),
    median=NA,
    mode=c(rep(p$umsy, nrow(df)),
          rep(pac$umsy, nrow(df)),
           ptva$umsy,
          rep(ptvb$umsy, nrow(df)),
          ptvab$umsy,
          phmma$umsy[phmma$regime],
          rep(phmmb$umsy,nrow(df)),
          phmm$umsy[phmm$regime]),
     convergence=rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence+ ptva$conv_problem,
      ptvb$model$convergence+ ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem),each=nrow(df)))

    dfumsy$pbias<- ((dfumsy$mode-dfumsy$sim)/dfumsy$sim)*100
    dfumsy$bias<- (dfumsy$mode-dfumsy$sim)


     #AIC
    dfaic<- data.frame(parameter="AIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",8),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma", "hmmb","hmmab"),
                       by=rep(NA,8),
                       sim=rep(NA,8),
                       median=NA,
                       mode=c(p$AICc,
                             pac$AICc,
                             ptva$AICc,
                             ptvb$AICc,
                             ptvab$AICc,
                             phmma$AICc,
                             phmmb$AICc,
                             phmm$AICc),
                       convergence=c(p$model$convergence + p$conv_problem,
                                     pac$model$convergence + pac$conv_problem,
                                     ptva$model$convergence+ ptva$conv_problem,
                                     ptvb$model$convergence+ ptvb$conv_problem,
                                     ptvab$model$convergence + ptvab$conv_problem,
                                     phmma$model$convergence + phmma$conv_problem,
                                     phmmb$model$convergence + phmmb$conv_problem,
                                     phmm$model$convergence + phmm$conv_problem),
                       pbias=rep(NA,8),
                       bias=rep(NA,8))
    #BIC
    dfbic<- data.frame(parameter="BIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",8),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma", "hmmb","hmmab"),
                       by=rep(NA,8),
                       sim=rep(NA,8),
                       median=NA,
                       mode=c(p$BIC,
                             pac$BIC,
                             ptva$BIC,
                             ptvb$BIC,
                             ptvab$BIC,
                             phmma$BIC,
                             phmmb$BIC,
                             phmm$BIC),
                       convergence=c(p$model$convergence + p$conv_problem,
                                     pac$model$convergence + pac$conv_problem,
                                     ptva$model$convergence+ ptva$conv_problem,
                                     ptvb$model$convergence+ ptvb$conv_problem,
                                     ptvab$model$convergence + ptvab$conv_problem,
                                     phmma$model$convergence + phmma$conv_problem,
                                     phmmb$model$convergence + phmmb$conv_problem,
                                     phmm$model$convergence + phmm$conv_problem),
                       pbias=rep(NA,8),
                       bias=rep(NA,8))
  
   #lfo
    lfostatic <- tmb_mod_lfo_cv(data=df,model='static', L=15)
    lfoac <- tmb_mod_lfo_cv(data=df,model='staticAC', L=15)
    lfoalpha <- tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=15)
    lfobeta <- tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=15)
    lfoalphabeta <- tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=15)
    lfohmma <- tmb_mod_lfo_cv(data=df,model='HMM_a', L=15, dirichlet_prior=dirpr)
    lfohmmb <- tmb_mod_lfo_cv(data=df,model='HMM_b', L=15, dirichlet_prior=dirpr)
    lfohmm <- tmb_mod_lfo_cv(data=df,model='HMM', L=15, dirichlet_prior=dirpr)
    
    dflfo<- data.frame(parameter="LFO",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",20),
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
                        median=NA,
                       mode=c(sum(lfostatic$lastparam), 
                           sum(lfoac$lastparam), 
                           sum(lfoalpha$lastparam), 
                           sum(lfoalpha$last3paramavg), 
                           sum(lfoalpha$last5paramavg), 
                           sum(lfobeta$lastparam), 
                           sum(lfobeta$last3paramavg), 
                           sum(lfobeta$last5paramavg), 
                           sum(lfoalphabeta$lastparam), 
                           sum(lfoalphabeta$last3paramavg), 
                           sum(lfoalphabeta$last5paramavg),    
                           sum(lfohmma$lastregime_pick),
                           sum(lfohmma$last3regime_pick),
                           sum(lfohmma$last5regime_pick),
                           sum(lfohmmb$lastregime_pick),
                           sum(lfohmmb$last3regime_pick),
                           sum(lfohmmb$last5regime_pick),
                           sum(lfohmm$lastregime_pick),
                           sum(lfohmm$last3regime_pick),
                           sum(lfohmm$last5regime_pick)
                           ),
                       convergence=c(sum(lfostatic$conv_problem),
                           sum(lfoac$conv_problem), 
                           sum(lfoalpha$conv_problem), 
                           sum(lfoalpha$conv_problem), 
                           sum(lfoalpha$conv_problem), 
                           sum(lfobeta$conv_problem), 
                           sum(lfobeta$conv_problem), 
                           sum(lfobeta$conv_problem), 
                           sum(lfoalphabeta$conv_problem), 
                           sum(lfoalphabeta$conv_problem), 
                           sum(lfoalphabeta$conv_problem),    
                           sum(lfohmma$conv_problem),
                           sum(lfohmma$conv_problem),
                           sum(lfohmma$conv_problem),
                           sum(lfohmmb$conv_problem),
                           sum(lfohmmb$conv_problem),
                           sum(lfohmmb$conv_problem),
                           sum(lfohmm$conv_problem),
                           sum(lfohmm$conv_problem),
                           sum(lfohmm$conv_problem)
                                     ),
                       pbias=rep(NA,20),
                       bias=rep(NA,20))



   
    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic,dflfo)

  return(dff)

}
