
tmb_func <- function(path=".",a, u) {
    
  allsimest <- list()
  simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))


  p <- ricker_TMB(data=df)
  pac <- ricker_TMB(data=df, AC=TRUE)
  ptva <- ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- ricker_rw_TMB(data=df, tv.par='both')
  phmma <- ricker_hmm_TMB(data=df, tv.par='a')
  phmmb <- ricker_hmm_TMB(data=df, tv.par='b')
  phmm  <- ricker_hmm_TMB(data=df, tv.par='both')


  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",8)),each=nrow(df)),
              model=rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
              by=rep(dat$year,8),
              sim=rep(dat$alpha,8),
              est=c(rep(p$alpha,nrow(df)),
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
                    
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",8)),each=nrow(df)),
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(1/dat$beta,8),
      est=c(rep(p$Smax,nrow(df)),
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
      
    dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100

       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",8)),each=nrow(df)),
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(dat$sigma,8),
      est=c(rep(p$sig,nrow(df)),
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
    dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100

              
    #Smsy
    smsysim<-smsyCalc(dat$alpha,dat$beta)
  
    dfsmsy<- data.frame(parameter="smsy",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",8)),each=nrow(df)),
      model=rep(c("simple",
        "autocorr",
        "rwa","rwb","rwab",
        "hmma_regime", "hmmb_regime", "hmmab_regime"),each=nrow(df)),
      by=rep(dat$year,8),
      sim=rep(smsysim,8),
      est=c(rep(p$Smsy,nrow(df)),
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
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",8)),each=nrow(df)),
    model=rep(c("simple",
      "autocorr",
      "rwa","rwb","rwab",
      "hmma_regime", "hmmb_regime", "hmmab_regime"),each=nrow(df)),
    by=rep(dat$year,8),
    sim=rep(unlist(mapply(sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),8),
    est=c(unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="simple"&dfa$method=="MLE"],
            Smsy=dfsmsy$est[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
            b=1/dfsmax$est[dfsmax$model=="simple"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="autocorr"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="autocorr"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwa"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="rwa"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwb"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="rwb"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwab"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="rwab"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmma_regime"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="hmma_regime"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"&dfsmsy$method=="MLE"],
           b=1/dfsmax$est[dfsmax$model=="hmmb_regime"&dfsmax$method=="MLE"])),
        unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"&dfa$method=="MLE"],
          Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$est[dfsmax$model=="hmmab_regime"&dfsmax$method=="MLE"]))),
     convergence=rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptvab$conv_problem,
      ptvb$model$convergence + ptvab$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem),
      each=nrow(df)))
  
    dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100

         
  #umsy
  
  dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",8)),each=nrow(df)),
    model=rep(c("simple",
      "autocorr",
      "rwa","rwb","rwab",
      "hmma_regime", "hmmb_regime","hmmab_regime"),each=nrow(df)),
    by=rep(dat$year,8),
    sim=rep(umsyCalc(dat$alpha),8),
    est=c(rep(p$umsy, nrow(df)),
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

    dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100

    #AIC
    dfaic<- data.frame(parameter="AIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",8),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime", "hmmb_regime","hmmab_regime"),
                       by=rep(NA,8),
                       sim=rep(NA,8),
                       est=c(p$AICc,
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
                       pbias=rep(NA,8))

    #BIC
    dfbic<- data.frame(parameter="BIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",8),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime", "hmmb_regime","hmmab_regime"),
                       by=rep(NA,8),
                       sim=rep(NA,8),
                       est=c(p$BIC,
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
                       pbias=rep(NA,8))
    
    #new pointwise LL for AIC variants
    LL_mat=matrix(ncol=nrow(df),nrow=5)
    
    LL_mat[1,]=log(dnorm(df$logRS,mean=p$alpha-p$beta*df$S,sd=p$sig))
    LL_mat[2,1]=log(dnorm(df$logRS[1],mean=pac$alpha-pac$beta*df$S[1],sd=pac$sig))
    LL_mat[2,2:nrow(df)]=log(dnorm(df$logRS[-1],mean=pac$alpha-pac$beta*df$S[-1]+pac$rho*pac$residuals[-nrow(df)],sd=pac$sigar))
    LL_mat[3,]=log(dnorm(df$logRS,mean=ptva$alpha-ptva$beta*df$S,sd=ptva$sig))
    LL_mat[4,]=log(dnorm(df$logRS,mean=ptvb$alpha-ptvb$beta*df$S,sd=ptvb$sig))
    LL_mat[5,]=log(dnorm(df$logRS,mean=ptvab$alpha-ptvab$beta*df$S,sd=ptvab$sig))
 
    LL_matd90=LL_mat[,apply(LL_mat,2,sum)>=quantile(apply(LL_mat,2,sum),0.1)]
    LL_matd80=LL_mat[,apply(LL_mat,2,sum)>=quantile(apply(LL_mat,2,sum),0.2)]
    
    npar=c(3,4,4,4,5)
   
    
    nll=-apply(LL_mat,1,sum)
    nlld90=-apply(LL_matd90,1,sum)
    nlld80=-apply(LL_matd80,1,sum)
    
    #normal AIC to compare
    aic_n=2*nll + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
    bic_n= 2*nll + npar*log(nrow(df))
    
    dfaic2<- data.frame(parameter="AIC_n",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",5),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab"),
                       by=rep(NA,5),
                       sim=rep(NA,5),
                       est=aic_n,
                       convergence=c(p$model$convergence + p$conv_problem,
                                     pac$model$convergence + pac$conv_problem,
                                     ptva$model$convergence+ ptva$conv_problem,
                                     ptvb$model$convergence+ ptvb$conv_problem,
                                     ptvab$model$convergence + ptvab$conv_problem),
                       pbias=rep(NA,5))

     dfbic2<- data.frame(parameter="bic_n",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=bic_n,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
    

    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic,dfaic2,dfbic2)

  return(dff)

}



tmb_lfo_func <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  
  lfostatic<-samEst::tmb_mod_lfo_cv(data=df,model='static', L=10)
  lfoac <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='staticAC', L=10),error = function(e) {lfoac=list(lastparam=rep(-999,length(lfoac$lastparam)))})
  lfoalpha <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=10),error = function(e) {lfoalpha=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last5param=rep(-999,length(lfoac$lastparam)))})
  lfobeta <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=10),error = function(e) {lfobeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last5param=rep(-999,length(lfoac$lastparam)))})
  lfoalphabeta <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=10),error = function(e) {lfoalphabeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last5param=rep(-999,length(lfoac$lastparam)))})
  LLdf<-rbind(lfostatic$lastparam,lfoac$lastparam,
              lfoalpha$lastparam,lfoalpha$last3param,lfoalpha$last5param,
              lfobeta$lastparam,lfobeta$last3param,lfobeta$last5param,
              lfoalphabeta$lastparam,lfoalphabeta$last3param,lfoalphabeta$last5param
  )
 
  return(LLdf)
  
}