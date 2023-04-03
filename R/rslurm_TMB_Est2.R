
###Tested and working TMB pbias & AIC/BIC"####
library(rslurm)

library(samEst,lib="/fs/vnas_Hdfo/comda/dag004/Rlib")


simPars <- read.csv("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/data/generic/SimPars.csv")
#save(simPars, file = "data/harcnkSimPars.RData")
#load("data/harcnkSimPars.RData")


tmb_func <- function(a,u) {
    
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  p <- samEst::ricker_TMB(data=df)
  pac <- samEst::ricker_TMB(data=df, AC=TRUE)
  ptva <- samEst::ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- samEst::ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- samEst::ricker_rw_TMB(data=df, tv.par='both')
  phmma <- samEst::ricker_hmm_TMB(data=df, tv.par='a')
  phmmb <- samEst::ricker_hmm_TMB(data=df, tv.par='b')
  phmm  <- samEst::ricker_hmm_TMB(data=df, tv.par='both')


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
    smsysim<-samEst::smsyCalc(dat$alpha,dat$beta)
  
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
                       sim=rep(unlist(mapply(samEst::sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),8),
                       est=c(unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="simple"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="simple"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="autocorr"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="autocorr"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwa"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwa"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwb"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwb"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwab"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwab"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmma_regime"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmma_regime"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"&dfsmsy$method=="MLE"],
                                           b=1/dfsmax$est[dfsmax$model=="hmmb_regime"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"&dfa$method=="MLE"],
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
    npar2=c(3,4,3+log(nrow(df)),3+log(nrow(df)),3+2*log(nrow(df)))
    
    nll=-apply(LL_mat,1,sum)
    nlld90=-apply(LL_matd90,1,sum)
    nlld80=-apply(LL_matd80,1,sum)
    
    #normal AIC to compare
    aic_n=2*nll + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
    bic_n= 2*nll + npar*log(nrow(df))
    #d90
    aic_d90=2*nlld90 + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
    bic_d90= 2*nlld90 + npar*log(nrow(df))
    #d80
    aic_d80=2*nlld80 + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
    bic_d80= 2*nlld80 + npar*log(nrow(df))
    
    #with revised number of parameters
    aic_n2=2*nll + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
    bic_n2= 2*nll + npar2*log(nrow(df))
    #d90
    aic_d902=2*nlld90 + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
    bic_d902= 2*nlld90 + npar2*log(nrow(df))
    #d80
    aic_d802=2*nlld80 + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
    bic_d802= 2*nlld80 + npar2*log(nrow(df))
    
    
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
    
    dfaicd90<- data.frame(parameter="AIC_d90",
                          iteration=u,
                          scenario= simPars$scenario[a],
                          method=rep("MLE",5),
                          model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),
                          by=rep(NA,5),
                          sim=rep(NA,5),
                          est=aic_d90,
                          convergence=c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),
                          pbias=rep(NA,5))
    
    dfaicd80<-data.frame(parameter="AIC_d80",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=aic_d80,
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
    
    dfbicd90<- data.frame(parameter="bic_d90",
                          iteration=u,
                          scenario= simPars$scenario[a],
                          method=rep("MLE",5),
                          model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),
                          by=rep(NA,5),
                          sim=rep(NA,5),
                          est=bic_d90,
                          convergence=c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),
                          pbias=rep(NA,5))
    
    dfbicd80<-data.frame(parameter="bic_d80",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=bic_d80,
                         convergence=c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence+ ptva$conv_problem,
                                       ptvb$model$convergence+ ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),
                         pbias=rep(NA,5))
    
    dfaicnpar2<- data.frame(parameter="AIC_npar2",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=aic_n2,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
    dfaicd902<- data.frame(parameter="AIC_d90npar2",
                          iteration=u,
                          scenario= simPars$scenario[a],
                          method=rep("MLE",5),
                          model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),
                          by=rep(NA,5),
                          sim=rep(NA,5),
                          est=aic_d902,
                          convergence=c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),
                          pbias=rep(NA,5))
    
    dfaicd802<-data.frame(parameter="AIC_d80npar2",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=aic_d80,
                         convergence=c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence+ ptva$conv_problem,
                                       ptvb$model$convergence+ ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),
                         pbias=rep(NA,5))
    
    dfbic22<- data.frame(parameter="bic_npar2",
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
    
    dfbicd902<- data.frame(parameter="bic_d90npar2",
                          iteration=u,
                          scenario= simPars$scenario[a],
                          method=rep("MLE",5),
                          model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),
                          by=rep(NA,5),
                          sim=rep(NA,5),
                          est=bic_d90,
                          convergence=c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),
                          pbias=rep(NA,5))
    
    dfbicd802<-data.frame(parameter="bic_d80npar2",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=bic_d80,
                         convergence=c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence+ ptva$conv_problem,
                                       ptvb$model$convergence+ ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),
                         pbias=rep(NA,5))
    
    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic,dfaic2,dfaicd90,dfaicd80,dfbic2,dfbicd90,dfbicd80,dfaicnpar2,dfaicd902,dfaicd802,dfbic22,dfbicd902,dfbicd802)
    
  return(dff)

}


tmb_lfo_func <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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

pars<-data.frame(a=rep(1:12,each=2000),
  u=1:2000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 1, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                    libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                    global_objects=c("simPars"))

save.image("./slrmjb.RData")
pars<-data.frame(a=rep(1:12,each=5),
                 u=1:5)
#load again
library(rslurm)
setwd('..')
load('slrmjb.RData')

res <- get_slurm_out(sjobtmb, outtype = 'table', wait = FALSE)
head(res, 3)

sc14=subset(res, scenario %in% c('stationary','autocorr','decLinearProd','regimeProd'))
sc58=subset(res, scenario %in% c('sineProd','regimeCap','decLinearCap','sigmaShift'))
sc912=subset(res, scenario %in% c('regimeProdCap','shiftCap','decLinearProdshiftCap','shiftProd'))

saveRDS(sc14,file='sc1_4.RDS')
saveRDS(sc58,file='sc5_8.RDS')
saveRDS(sc912,file='sc9_12.RDS')
#Stan
cmdstanr::set_cmdstan_path(path="/fs/vnas_Hdfo/comda/dag004/.cmdstan/cmdstan-2.31.0")
file1=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1f.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2f.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f.stan")
mod3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f.stan")
mod5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6f.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7f.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8f.stan")
mod8=cmdstanr::cmdstan_model(file8)


stan_func<- function(a,u){
  cmdstanr::set_cmdstan_path(path="/fs/vnas_Hdfo/comda/dag004/.cmdstan/cmdstan-2.31.0")
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
             alpha_dirichlet=c(1,1)
  )
  
  #create folder to hold temp files
#  dir.create(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u,sep=''))
  
#  ls=list.files("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/")
#  if(length(ls)>5){
 #   unlink(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u-2,
#                 '/*',sep=''))
#}
  
  #
  f1 <- mod1$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f2 <- mod2$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f4 <- mod4$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f5 <- mod5$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f6 <- mod6$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f7 <- mod7$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f8 <- mod8$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  #Max. prod
  
  #regime state sequence:
  phmma_alpha=f6$summary(variables=c('log_a'),'median')$median
  phmma_alpha_regime=phmma_alpha[f6$summary(variables=c('zstar'),'median')$median]
  phmmab_alpha=f8$summary(variables=c('log_a'),'median')$median
  phmmab_alpha_regime=phmmab_alpha[f8$summary(variables=c('zstar'),'median')$median]
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   est=c(rep(f1$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         rep(f2$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f3$summary(variables=c('log_a'),'median')$median,
                         rep(f4$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f5$summary(variables=c('log_a'),'median')$median,
                         phmma_alpha_regime,
                         rep(f7$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         phmmab_alpha_regime
                   ))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
 
  #Smax
  phmmb_smax=f7$summary(variables=c('S_max'),'median')$median
  phmmb_smax_regime=phmmb_smax[f7$summary(variables=c('zstar'),'median')$median]
  phmmab_smax=f8$summary(variables=c('S_max'),'median')$median
  phmmab_smax_regime=phmmb_smax[f8$summary(variables=c('zstar'),'median')$median]
  
  dfsmax<- data.frame(parameter="smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(1/dat$beta,8),
                      est=c(rep(f1$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f2$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f3$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            f4$summary(variables=c('S_max'),'median')$median,
                            f5$summary(variables=c('S_max'),'median')$median,
                            rep(f6$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            phmmb_smax_regime,
                            phmmab_smax_regime
                      ))
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  #obs error
  dfsig<- data.frame(parameter="sigma_obs",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma_regime","hmmb_regime","hmmab_regime"),
                     by=NA,
                     sim=NA,
                     est=c(f1$summary(variables=c('sigma'),'median')$median,
                           f2$summary(variables=c('sigma'),'median')$median,
                           f3$summary(variables=c('sigma'),'median')$median,
                           f4$summary(variables=c('sigma'),'median')$median,
                           f5$summary(variables=c('sigma'),'median')$median,
                           f6$summary(variables=c('sigma'),'median')$median,
                           f7$summary(variables=c('sigma'),'median')$median,
                           f8$summary(variables=c('sigma'),'median')$median))
  
  dfsig$pbias<- NA
  
  #sigma a
  dfsiga<- data.frame(parameter="sigma_a",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwa","rwab"),
                      by=NA,
                      sim=NA,
                      est=c(f3$summary(variables=c('sigma_a'),'median')$median,
                            f5$summary(variables=c('sigma_a'),'median')$median))
  
  dfsiga$pbias<- NA
  
  #sigma b
  dfsigb<- data.frame(parameter="sigma_b",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwb","rwab"),
                      by=NA,
                      sim=NA,
                      est=c(f4$summary(variables=c('sigma_b'),'median')$median,
                            f5$summary(variables=c('sigma_b'),'median')$median))
  
  dfsigb$pbias<- NA
  #S msy
  #Smsy - estimate from posterior 
  #static
  smsysim<-samEst::smsyCalc(dat$alpha,dat$beta)
  
  dfsmsy<- data.frame(parameter="smsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(smsysim,8),
                      est=c(samEst::smsyCalc(a=dfa$est[dfa$model=="simple"],b=1/dfsmax$est[dfsmax$model=="simple"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="autocorr"],b=1/dfsmax$est[dfsmax$model=="autocorr"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwa"],b=1/dfsmax$est[dfsmax$model=="rwa"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwb"],b=1/dfsmax$est[dfsmax$model=="rwb"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwab"],b=1/dfsmax$est[dfsmax$model=="rwab"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmma_regime"],b=1/dfsmax$est[dfsmax$model=="hmma_regime"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmb_regime"],b=1/dfsmax$est[dfsmax$model=="hmmb_regime"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmab_regime"],b=1/dfsmax$est[dfsmax$model=="hmmab_regime"])))
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100
  
  
  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                       model=rep(c("simple",
                                   "autocorr",
                                   "rwa","rwb","rwab",
                                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                       by=rep(dat$year,8),
                       sim=rep(unlist(mapply(samEst::sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),8),
                       est=c(unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="simple"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="simple"], 
                                           b=1/dfsmax$est[dfsmax$model=="simple"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="autocorr"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="autocorr"], 
                                           b=1/dfsmax$est[dfsmax$model=="autocorr"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwa"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwa"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwa"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwb"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwb"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwab"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwab"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwab"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmma_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmma_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmb_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmab_regime"]))))
  
  dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100     
  
  #umsy
  dfumsy<- data.frame(parameter="umsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(samEst::umsyCalc(dat$alpha),8),
                      est=c(samEst::umsyCalc(dfa$est[dfa$model=="simple"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="autocorr"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwa"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwb"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwab"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmma_regime"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmb_regime"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmab_regime"]))
  )
  
  dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100
  
  #Pointwise loglikelihoods
  dfelpd<- data.frame(parameter="ELPD",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",8),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(NA,8),
                      est=c(
                        f1$summary(variables=c('log_lik'),'median')$median,
                        f2$summary(variables=c('log_lik'),'median')$median,
                        f3$summary(variables=c('log_lik'),'median')$median,
                        f4$summary(variables=c('log_lik'),'median')$median,
                        f5$summary(variables=c('log_lik'),'median')$median,
                        f6$summary(variables=c('log_lik'),'median')$median,
                        f7$summary(variables=c('log_lik'),'median')$median,
                        f8$summary(variables=c('log_lik'),'median')$median)
                      ,
                      pbias=rep(NA,8))
  
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsiga,dfsigb,dfsmsy,dfsgen,dfumsy,dfelpd)
  
  return(dff)
  
}

pars1<-data.frame(a=rep(1:2,each=100),
                  u=1:100)

pars2<-data.frame(a=rep(5:8,each=1000),
                 u=1:1000)

pars3<-data.frame(a=rep(9:12,each=1000),
                  u=1:1000)

sjobstan1 <- slurm_apply(stan_func, pars1, jobname = 'Srun1',
                         nodes = 1, cpus_per_node = 12, submit = FALSE,
                         pkgs=c("samEst","cmdstanr"),
                         rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                         libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                         global_objects=c("simPars","mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"))

sjobstan2 <- slurm_apply(stan_func, pars2, jobname = 'Srun2',
                         nodes = 1, cpus_per_node = 12, submit = FALSE,
                         pkgs=c("samEst","cmdstanr"),
                         rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                         libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                         global_objects=c("simPars","mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"))

sjobstan3 <- slurm_apply(stan_func, pars3, jobname = 'Srun3',
                         nodes = 1, cpus_per_node = 12, submit = FALSE,
                         pkgs=c("samEst","cmdstanr"),
                         rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                         libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                         global_objects=c("simPars","mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"))

save.image("./slrmjb.RData")

q()


library(rslurm)
setwd('..')
load('slrmjb.RData')

res1 <- get_slurm_out(sjobstan1, outtype = 'table', wait = FALSE)
head(res1, 3)

res <- get_slurm_out(sjobstan2, outtype = 'table', wait = FALSE)
head(res, 3)

res <- get_slurm_out(sjobstan3, outtype = 'table', wait = FALSE)
head(res, 3)


#https://portal.science.gc.ca/confluence/display/SCIDOCS/Quick+Start+to+Using+Linux+Clusters+With+SLURM



#Productivity scenarios####
library(rslurm)

library(samEst,lib="/fs/vnas_Hdfo/comda/dag004/Rlib")

simPars <- data.frame(scenario=c('trendLinearProd1','trendLinearProd2','trendLinearProd5','trendLinearProd7','regimeProd1','regimeProd2','regimeProd5','regimeProd7'),nameOM=c('trendLinearProd1','trendLinearProd2','trendLinearProd5','trendLinearProd7','regimeProd1','regimeProd2','regimeProd5','regimeProd7'),nameMP='fixedER')

tmb_func_prod <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/prod_scns/",simPars$scenario[a],"/",simPars$nameOM[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                   S=dat$obsSpawners,
                   R=dat$obsRecruits,
                   logRS=log(dat$obsRecruits/dat$obsSpawners))
  
  p <- samEst::ricker_TMB(data=df)
  pac <- samEst::ricker_TMB(data=df, AC=TRUE)
  ptva <- samEst::ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- samEst::ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- samEst::ricker_rw_TMB(data=df, tv.par='both')
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MLE",5)),each=nrow(df)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa",'rwb','rwab'),each=nrow(df)),
                   by=rep(dat$year,5),
                   sim=rep(dat$alpha,5),
                   est=c(rep(p$alpha,nrow(df)),
                         rep(pac$alpha,nrow(df)),
                         ptva$alpha,
                         rep(ptvb$alpha,nrow(df)),
                         ptvab$alpha),
                   convergence=c(rep(c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence + ptva$conv_problem,
                                       ptvb$model$convergence + ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem
                   ),each=nrow(df))))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MLE",5)),each=nrow(df)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),each=nrow(df)),
                      by=rep(dat$year,5),
                      sim=rep(1/dat$beta,5),
                      est=c(rep(p$Smax,nrow(df)),
                            rep(pac$Smax,nrow(df)),
                            rep(ptva$Smax,nrow(df)),
                            ptvb$Smax,
                            ptvab$Smax),
                      convergence=rep(c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence + ptva$conv_problem,
                                        ptvb$model$convergence + ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),each=nrow(df)))
  
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  
  #sigma
  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep(c(rep("MLE",5)),each=nrow(df)),
                     model=rep(c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),each=nrow(df)),
                     by=rep(dat$year,5),
                     sim=rep(dat$sigma,5),
                     est=c(rep(p$sig,nrow(df)),
                           rep(pac$sig,nrow(df)),
                           rep(ptva$sig,nrow(df)),
                           rep(ptvb$sig,nrow(df)),
                           rep(ptvab$sig,nrow(df))
                           ),
                     convergence=rep(c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence + ptva$conv_problem,
                                       ptvb$model$convergence + ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),each=nrow(df)))
  dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
  
  #Smsy
  smsysim<-samEst::smsyCalc(dat$alpha,dat$beta)
  
  dfsmsy<- data.frame(parameter="smsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MLE",5)),each=nrow(df)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),each=nrow(df)),
                      by=rep(dat$year,5),
                      sim=rep(smsysim,5),
                      est=c(rep(p$Smsy,nrow(df)),
                            rep(pac$Smsy,nrow(df)),
                            ptva$Smsy,
                            ptvb$Smsy,
                            ptvab$Smsy),    
                      convergence=rep(c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence + ptva$conv_problem,
                                        ptvb$model$convergence + ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),each=nrow(df))) 
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100
  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MLE",5)),each=nrow(df)),
                       model=rep(c("simple",
                                   "autocorr",
                                   "rwa","rwb","rwab"),each=nrow(df)),
                       by=rep(dat$year,5),
                       sim=rep(unlist(mapply(samEst::sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),5),
                       est=c(unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="simple"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="simple"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="autocorr"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="autocorr"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwa"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwa"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwb"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwb"&dfsmax$method=="MLE"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwab"&dfa$method=="MLE"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwab"&dfsmax$method=="MLE"]))),
                       convergence=rep(c(p$model$convergence + p$conv_problem,
                                         pac$model$convergence + pac$conv_problem,
                                         ptva$model$convergence + ptvab$conv_problem,
                                         ptvb$model$convergence + ptvab$conv_problem,
                                         ptvab$model$convergence + ptvab$conv_problem),
                                       each=nrow(df)))
  
  dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100     
  #umsy
  
  dfumsy<- data.frame(parameter="umsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MLE",5)),each=nrow(df)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),each=nrow(df)),
                      by=rep(dat$year,5),
                      sim=rep(umsyCalc(dat$alpha),5),
                      est=c(rep(p$umsy, nrow(df)),
                            rep(pac$umsy, nrow(df)),
                            ptva$umsy,
                            rep(ptvb$umsy, nrow(df)),
                            ptvab$umsy),
                      convergence=rep(c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),each=nrow(df)))
  
  dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100
  
  #AIC
  dfaic<- data.frame(parameter="AIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MLE",5),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab"),
                     by=rep(NA,5),
                     sim=rep(NA,5),
                     est=c(p$AICc,
                           pac$AICc,
                           ptva$AICc,
                           ptvb$AICc,
                           ptvab$AICc),
                     convergence=c(p$model$convergence + p$conv_problem,
                                   pac$model$convergence + pac$conv_problem,
                                   ptva$model$convergence+ ptva$conv_problem,
                                   ptvb$model$convergence+ ptvb$conv_problem,
                                   ptvab$model$convergence + ptvab$conv_problem),
                     pbias=rep(NA,5))
  #BIC
  dfbic<- data.frame(parameter="BIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MLE",5),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab"),
                     by=rep(NA,5),
                     sim=rep(NA,5),
                     est=c(p$BIC,
                           pac$BIC,
                           ptva$BIC,
                           ptvb$BIC,
                           ptvab$BIC),
                     convergence=c(p$model$convergence + p$conv_problem,
                                   pac$model$convergence + pac$conv_problem,
                                   ptva$model$convergence+ ptva$conv_problem,
                                   ptvb$model$convergence+ ptvb$conv_problem,
                                   ptvab$model$convergence + ptvab$conv_problem),
                     pbias=rep(NA,5))
  
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
  npar2=c(3,4,3+log(nrow(df)),3+log(nrow(df)),3+2*log(nrow(df)))
  
  nll=-apply(LL_mat,1,sum)
  nlld90=-apply(LL_matd90,1,sum)
  nlld80=-apply(LL_matd80,1,sum)
  
  #normal AIC to compare
  aic_n=2*nll + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
  bic_n= 2*nll + npar*log(nrow(df))
  #d90
  aic_d90=2*nlld90 + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
  bic_d90= 2*nlld90 + npar*log(nrow(df))
  #d80
  aic_d80=2*nlld80 + 2*npar +(2*npar*(npar+1)/(nrow(df)-npar-1))
  bic_d80= 2*nlld80 + npar*log(nrow(df))
  
  #with revised number of parameters
  aic_n2=2*nll + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
  bic_n2= 2*nll + npar2*log(nrow(df))
  #d90
  aic_d902=2*nlld90 + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
  bic_d902= 2*nlld90 + npar2*log(nrow(df))
  #d80
  aic_d802=2*nlld80 + 2*npar2 +(2*npar2*(npar2+1)/(nrow(df)-npar2-1))
  bic_d802= 2*nlld80 + npar2*log(nrow(df))
  
  
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
  
  dfaicd90<- data.frame(parameter="AIC_d90",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=aic_d90,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
  
  dfaicd80<-data.frame(parameter="AIC_d80",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",5),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab"),
                       by=rep(NA,5),
                       sim=rep(NA,5),
                       est=aic_d80,
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
  
  dfbicd90<- data.frame(parameter="bic_d90",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=bic_d90,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
  
  dfbicd80<-data.frame(parameter="bic_d80",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",5),
                       model=c("simple",
                               "autocorr",
                               "rwa","rwb","rwab"),
                       by=rep(NA,5),
                       sim=rep(NA,5),
                       est=bic_d80,
                       convergence=c(p$model$convergence + p$conv_problem,
                                     pac$model$convergence + pac$conv_problem,
                                     ptva$model$convergence+ ptva$conv_problem,
                                     ptvb$model$convergence+ ptvb$conv_problem,
                                     ptvab$model$convergence + ptvab$conv_problem),
                       pbias=rep(NA,5))
  
  dfaicnpar2<- data.frame(parameter="AIC_npar2",
                          iteration=u,
                          scenario= simPars$scenario[a],
                          method=rep("MLE",5),
                          model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab"),
                          by=rep(NA,5),
                          sim=rep(NA,5),
                          est=aic_n2,
                          convergence=c(p$model$convergence + p$conv_problem,
                                        pac$model$convergence + pac$conv_problem,
                                        ptva$model$convergence+ ptva$conv_problem,
                                        ptvb$model$convergence+ ptvb$conv_problem,
                                        ptvab$model$convergence + ptvab$conv_problem),
                          pbias=rep(NA,5))
  dfaicd902<- data.frame(parameter="AIC_d90npar2",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=aic_d902,
                         convergence=c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence+ ptva$conv_problem,
                                       ptvb$model$convergence+ ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),
                         pbias=rep(NA,5))
  
  dfaicd802<-data.frame(parameter="AIC_d80npar2",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=aic_d80,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
  
  dfbic22<- data.frame(parameter="bic_npar2",
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
  
  dfbicd902<- data.frame(parameter="bic_d90npar2",
                         iteration=u,
                         scenario= simPars$scenario[a],
                         method=rep("MLE",5),
                         model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab"),
                         by=rep(NA,5),
                         sim=rep(NA,5),
                         est=bic_d90,
                         convergence=c(p$model$convergence + p$conv_problem,
                                       pac$model$convergence + pac$conv_problem,
                                       ptva$model$convergence+ ptva$conv_problem,
                                       ptvb$model$convergence+ ptvb$conv_problem,
                                       ptvab$model$convergence + ptvab$conv_problem),
                         pbias=rep(NA,5))
  
  dfbicd802<-data.frame(parameter="bic_d80npar2",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MLE",5),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab"),
                        by=rep(NA,5),
                        sim=rep(NA,5),
                        est=bic_d80,
                        convergence=c(p$model$convergence + p$conv_problem,
                                      pac$model$convergence + pac$conv_problem,
                                      ptva$model$convergence+ ptva$conv_problem,
                                      ptvb$model$convergence+ ptvb$conv_problem,
                                      ptvab$model$convergence + ptvab$conv_problem),
                        pbias=rep(NA,5))
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic,dfaic2,dfaicd90,dfaicd80,dfbic2,dfbicd90,dfbicd80,dfaicnpar2,dfaicd902,dfaicd802,dfbic22,dfbicd902,dfbicd802)
  
  return(dff)
  
}

pars<-data.frame(a=rep(1:8,each=10),
                 u=1:10)


sjobtmbprod <- slurm_apply(tmb_func, pars, jobname = 'TMBprodscn',
                       nodes = 1, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst"),
                       rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                       libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                       global_objects=c("simPars"))

save.image("./sjobtmbprod.RData")

#load again
library(rslurm)
setwd('..')
load('sjobtmbprod.RData')

res <- get_slurm_out(sjobtmb, outtype = 'table', wait = FALSE)
head(res, 3)







###LOCAL pbias & AIC/BIC"####
#test it out
tmb_func2 <- function(a, u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("./outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                   S=dat$obsSpawners,
                   R=dat$obsRecruits,
                   logRS=log(dat$obsRecruits/dat$obsSpawners))
  
  p <- ricker_TMB(data=df)
  pac <- samEst::ricker_TMB(data=df, AC=TRUE)
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
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic)
  
  return(dff)
  
}



simPars <- data.frame(scenario=c('trendLinearProd1','trendLinearProd2','trendLinearProd5','trendLinearProd7','regimeProd1','regimeProd2','regimeProd5','regimeProd7'),nameOM=c('trendLinearProd1','trendLinearProd2','trendLinearProd5','trendLinearProd7','regimeProd1','regimeProd2','regimeProd5','regimeProd7'),nameMP='fixedER')


df=list()
for(a in 1:8){
  df[[a]]=tmb_func2(a=a,u=1)
  for(u in 2:10){
    set=tmb_func2(a=a,u=u)
    df[[a]]=rbind(df[[a]],set)
    
  }
}

#summarize model selection
aic_set=list()
bic_set=list()
for(i in 1:length(df)){
  aic=subset(df[[i]],parameter=='AIC')
  aic_set[[i]]=tidyr::spread(aic[,-9],key=model,value=est)
  aic_set[[i]]=aic_set[[i]][c(15,8,12,14,13,9,11,10)]
  bic=subset(df[[i]],parameter=='BIC')
  bic_set[[i]]=tidyr::spread(bic[,-9],key=model,value=est)
  bic_set[[i]]=bic_set[[i]][c(15,8,12,14,13,9,11,10)]
}
sc1=apply(aic_set[[1]],1,which.min)
cn1=summary(factor(sc1),levels=seq(1:8))/1000
sc2=apply(aic_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:8)))/1000
sc3=apply(aic_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:8)))/1000
sc4=apply(aic_set[[7]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:8)))/1000
sc5=apply(aic_set[[11]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:8)))/1000
sc6=apply(aic_set[[4]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:8)))/1000
sc7=apply(aic_set[[6]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:8)))/1000
sc8=apply(aic_set[[9]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:8)))/1000

sx1=apply(bic_set[[1]],1,which.min)
ck1=summary(factor(sx1),levels=seq(1:8))/1000
sx2=apply(bic_set[[2]],1,which.min)
ck2=summary(factor(sx2,levels=seq(1:8)))/1000
sx3=apply(bic_set[[3]],1,which.min)
ck3=summary(factor(sx3,levels=seq(1:8)))/1000
sx4=apply(bic_set[[7]],1,which.min)
ck4=summary(factor(sx4,levels=seq(1:8)))/1000
sx5=apply(bic_set[[11]],1,which.min)
ck5=summary(factor(sx5,levels=seq(1:8)))/1000
sx6=apply(bic_set[[4]],1,which.min)
ck6=summary(factor(sx6,levels=seq(1:8)))/1000
sx7=apply(bic_set[[6]],1,which.min)
ck7=summary(factor(sx7,levels=seq(1:8)))/1000
sx8=apply(bic_set[[9]],1,which.min)
ck8=summary(factor(sx8,levels=seq(1:8)))/1000

##Confusion matrices
conf_matrix <-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab",
                              "regime.a", "regime.b","regime.ab"),OM=c("stat.",
                                                                       "ac.",
                                                                       "d.a","d.b","d.ab",
"r.a", "r.b","r.ab"))
                          
conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:8]=cn1
conf_matrix$w_AIC[9:16]=cn2
conf_matrix$w_AIC[17:24]=cn3
conf_matrix$w_AIC[25:32]=cn4
conf_matrix$w_AIC[33:40]=cn5
conf_matrix$w_AIC[41:48]=cn6
conf_matrix$w_AIC[49:56]=cn7
conf_matrix$w_AIC[57:64]=cn8
conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:8]=ck1
conf_matrix$w_BIC[9:16]=ck2
conf_matrix$w_BIC[17:24]=ck3
conf_matrix$w_BIC[25:32]=ck4
conf_matrix$w_BIC[33:40]=ck5
conf_matrix$w_BIC[41:48]=ck6
conf_matrix$w_BIC[49:56]=ck7
conf_matrix$w_BIC[57:64]=ck8

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

library(ggplot2)
p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")
p

b=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")

b

ps=cowplot::plot_grid(p,b,
          ncol=2,nrow=1,labels=c("A","B"),label_size=10)

ggsave(filename = "outs/conf_mat_tmb_aicbic.pdf",
       plot=ps,
       width=12,height=6)
ggsave(filename = "outs/conf_mat_tmb_aic.pdf",
       plot=p,
       width=8,height=6)

ggsave(filename = "outs/conf_mat_tmb_bic.pdf",
       plot=b,
       width=12,height=6)


#By model class
sc1m=ifelse(sc1<3,1,sc1)
sc1m=ifelse(sc1>2&sc1<6,2,sc1m)
sc1m=ifelse(sc1>5,3,sc1m)
cn1m=summary(factor(sc1m),levels=seq(1:3))/1000
sc2m=ifelse(sc2<3,1,sc2)
sc2m=ifelse(sc2>2&sc2<6,2,sc2m)
sc2m=ifelse(sc2>5,3,sc2m)
cn2m=summary(factor(sc2m),levels=seq(1:3))/1000
sc3m=ifelse(sc3<3,1,sc3)
sc3m=ifelse(sc3>2&sc3<6,2,sc3m)
sc3m=ifelse(sc3>5,3,sc3m)
cn3m=summary(factor(sc3m),levels=seq(1:3))/1000
sc4m=ifelse(sc4<3,1,sc4)
sc4m=ifelse(sc4>2&sc4<6,2,sc4m)
sc4m=ifelse(sc4>5,3,sc4m)
cn4m=summary(factor(sc4m),levels=seq(1:3))/1000
sc5m=ifelse(sc5<3,1,sc5)
sc5m=ifelse(sc5>2&sc5<6,2,sc5m)
sc5m=ifelse(sc5>5,3,sc5m)
cn5m=summary(factor(sc5m),levels=seq(1:3))/1000
sc6m=ifelse(sc6<3,1,sc6)
sc6m=ifelse(sc6>2&sc6<6,2,sc6m)
sc6m=ifelse(sc6>5,3,sc6m)
cn6m=summary(factor(sc6m),levels=seq(1:3))/1000
sc7m=ifelse(sc7<3,1,sc7)
sc7m=ifelse(sc7>2&sc7<6,2,sc7m)
sc7m=ifelse(sc7>5,3,sc7m)
cn7m=summary(factor(sc7m),levels=seq(1:3))/1000
sc8m=ifelse(sc8<3,1,sc8)
sc8m=ifelse(sc8>2&sc8<6,2,sc8m)
sc8m=ifelse(sc8>5,3,sc8m)
cn8m=summary(factor(sc8m),levels=seq(1:3))/1000

sx1m=ifelse(sx1<3,1,sx1)
sx1m=ifelse(sx1>2&sx1<6,2,sx1m)
sx1m=ifelse(sx1>5,3,sx1m)
ck1m=summary(factor(sx1m),levels=seq(1:3))/1000
sx2m=ifelse(sx2<3,1,sx2)
sx2m=ifelse(sx2>2&sx2<6,2,sx2m)
sx2m=ifelse(sx2>5,3,sx2m)
ck2m=summary(factor(sx2m),levels=seq(1:3))/1000
sx3m=ifelse(sx3<3,1,sx3)
sx3m=ifelse(sx3>2&sx3<6,2,sx3m)
sx3m=ifelse(sx3>5,3,sx3m)
ck3m=summary(factor(sx3m,levels=seq(1:3)))/1000
sx4m=ifelse(sx4<3,1,sx4)
sx4m=ifelse(sx4>2&sx4<6,2,sx4m)
sx4m=ifelse(sx4>5,3,sx4m)
ck4m=summary(factor(sx4m),levels=seq(1:3))/1000
sx5m=ifelse(sx5<3,1,sx5)
sx5m=ifelse(sx5>2&sx5<6,2,sx5m)
sx5m=ifelse(sx5>5,3,sx5m)
ck5m=summary(factor(sx5m),levels=seq(1:3))/1000
sx6m=ifelse(sx6<3,1,sx6)
sx6m=ifelse(sx6>2&sx6<6,2,sx6m)
sx6m=ifelse(sx6>5,3,sx6m)
ck6m=summary(factor(sx6m),levels=seq(1:3))/1000
sx7m=ifelse(sx7<3,1,sx7)
sx7m=ifelse(sx7>2&sx7<6,2,sx7m)
sx7m=ifelse(sx7>5,3,sx7m)
ck7m=summary(factor(sx7m),levels=seq(1:3))/1000
sx8m=ifelse(sx8<3,1,sx8)
sx8m=ifelse(sx8>2&sx8<6,2,sx8m)
sx8m=ifelse(sx8>5,3,sx8m)
ck8m=summary(factor(sx8m),levels=seq(1:3))/1000

##Confusion matrices
conf_matrix <-expand.grid(EM=c('static','dynamic','regime'),OM=c("stat",
                               "acorr",
                               "dyn.a","dyn.b","dyn.ab",
                               "reg.a", "reg.b","reg.ab"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:3]=cn1m
conf_matrix$w_AIC[4:6]=cn2m
conf_matrix$w_AIC[7:9]=cn3m
conf_matrix$w_AIC[10:12]=cn4m
conf_matrix$w_AIC[13:15]=cn5m
conf_matrix$w_AIC[16:18]=cn6m
conf_matrix$w_AIC[19:21]=cn7m
conf_matrix$w_AIC[22:24]=cn8m
conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:3]=ck1m
conf_matrix$w_BIC[4:6]=ck2m
conf_matrix$w_BIC[7:9]=ck3m
conf_matrix$w_BIC[10:12]=ck4m
conf_matrix$w_BIC[13:15]=ck5m
conf_matrix$w_BIC[16:18]=ck6m
conf_matrix$w_BIC[19:21]=ck7m
conf_matrix$w_BIC[22:24]=ck8m

library(ggplot2)
p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")
p

b=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")

b

ps=cowplot::plot_grid(p,b,
                      ncol=2,nrow=1,labels=c("A","B"),label_size=10)

ggsave(filename = "outs/conf_mat_mc_tmb_aicbic.pdf",
       plot=ps,
       width=12,height=4)
ggsave(filename = "outs/conf_mat_mc_tmb_aic.pdf",
       plot=p,
       width=8,height=5)

ggsave(filename = "outs/conf_mat_mc_tmb_bic.pdf",
       plot=b,
       width=12,height=5)



###Stan pbias & AIC/BIC####
#stan_func <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                   S=dat$obsSpawners,
                   R=dat$obsRecruits,
                   logRS=log(dat$obsRecruits/dat$obsSpawners))
  
  
  m1f=samEst::sr_mod2(type='static',ac = FALSE,par='n',lfo =F)
  m2f=samEst::sr_mod2(type='static',ac = TRUE,par='n',lfo=F)
  m3f=samEst::sr_mod2(type='rw',par='a',lfo=F)
  m4f=samEst::sr_mod2(type='rw',par='b',lfo=F)
  m5f=samEst::sr_mod2(type='rw',par='both',lfo=F)
  m6f=samEst::sr_mod2(type='hmm',par='a',lfo=F)
  m7f=samEst::sr_mod2(type='hmm',par='b',lfo=F)
  m8f=samEst::sr_mod2(type='hmm',par='both',lfo=F)
  #
  p <- samEst::ricker_stan(data=df,iter = 600, mod=m1f)
  #ricker autocorr
  pac <- samEst::ricker_stan(data=df,iter = 600, AC=TRUE, mod=m2f)
  #ricker tva
  ptva <- samEst::ricker_rw_stan(data=df, par="a",iter = 600, mod=m3f)
  #ricker tvb
  ptvb <- samEst::ricker_rw_stan(data=df, par="b",iter = 600, mod=m4f)
  #ricker tvab
  ptvab <- samEst::ricker_rw_stan(data=df, par="both",iter = 600, mod=m5f) 
  #ricker tvhmma
  phmma <- samEst::ricker_hmm_stan(data=df, par="a",iter = 600, mod=m6f)
  #ricker tvhmmb
  phmmb <- samEst::ricker_hmm_stan(data=df, par="b",iter = 600, mod=m7f)
  #ricker tvhmmab
  phmmab <- samEst::ricker_hmm_stan(data=df, par="both",iter = 600, mod=m8f) 
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",8)),each=nrow(df)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   est=c(rep(p$alpha,nrow(df)),
                         rep(pac$alpha,nrow(df)),
                         ptva$alpha[2:c(nrow(df)+1)],
                         rep(mean(btvb$alpha),nrow(df)),
                         ptva$alpha[2:c(nrow(df)+1)],
                         phmma$alpha_regime,
                         rep(phmmb$alpha,nrow(df)),
                         phmm$alpha_regime))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(df)),
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
                            phmm$Smax[phmm$regime]))
  
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  
  #sigma
  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep(c(rep("MCMC",8)),each=nrow(df)),
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
                           rep(phmm$sigma,nrow(df))))
  dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
  #AIC
  pll=rstan::extract(p)
  pacll=rstan::extract(pac)
  ptvall=rstan::extract(ptva)
  ptvbll=rstan::extract(ptvb)
  ptvabll=rstan::extract(ptvab)
  phmmall=rstan::extract(phmma)
  phmmbll=rstan::extract(phmmb)
  phmmabll=rstan::extract(phmm)
  
  dll=list(pll$log_lik,pacll$log_lik,ptvall$log_lik,ptvbll$log_lik,plltvab$log_lik,phmmall$log_lik,phmmbll$log_lik,phmmabll$log_lik)
  
  aic_weights=samEst::stan_aic(dll,form='aic',type='full',k=c(3,4,4,4,5,5,5,6))
  aic_weightsd90=samEst::stan_aic(dll,form='aic',type='d90',k=c(3,4,4,4,5,5,5,6))
  aic_weightsd80=samEst::stan_aic(dll,form='aic',type='full',k=c(3,4,4,4,5,5,5,6))
  bic_weights=samEst::stan_bic(dll,form='bic',type='full',k=c(3,4,4,4,5,5,5,6))
  bic_weightsd90=samEst::stan_bic(dll,form='bic',type='d90',k=c(3,4,4,4,5,5,5,6))
  bic_weightsd80=samEst::stan_bic(dll,form='bic',type='full',k=c(3,4,4,4,5,5,5,6))
  
  dfaic<- data.frame(parameter="AIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma_regime", "hmmb_regime","hmmab_regime"),
                     by=rep(NA,8),
                     sim=rep(NA,8),
                     est=c(aic_weights),
                     pbias=rep(NA,8))
  dfaicd90<- data.frame(parameter="AICd90",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma_regime", "hmmb_regime","hmmab_regime"),
                     by=rep(NA,8),
                     sim=rep(NA,8),
                     est=c(aic_weightsd90),
                     pbias=rep(NA,8))
  dfaicd80<- data.frame(parameter="AICd80",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MCMC",8),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab",
                                "hmma_regime", "hmmb_regime","hmmab_regime"),
                        by=rep(NA,8),
                        sim=rep(NA,8),
                        est=c(aic_weightsd80),
                        pbias=rep(NA,8))
  #BIC
  dfbic<- data.frame(parameter="BIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma_regime", "hmmb_regime","hmmab_regime"),
                     by=rep(NA,8),
                     sim=rep(NA,8),
                     est=c(bic_weights),
                     pbias=rep(NA,8))
  dfbicd90<- data.frame(parameter="BICd90",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MCMC",8),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab",
                                "hmma_regime", "hmmb_regime","hmmab_regime"),
                        by=rep(NA,8),
                        sim=rep(NA,8),
                        est=c(bic_weightsd90),
                        pbias=rep(NA,8))
  dfbicd80<- data.frame(parameter="BICd80",
                        iteration=u,
                        scenario= simPars$scenario[a],
                        method=rep("MCMC",8),
                        model=c("simple",
                                "autocorr",
                                "rwa","rwb","rwab",
                                "hmma_regime", "hmmb_regime","hmmab_regime"),
                        by=rep(NA,8),
                        sim=rep(NA,8),
                        est=c(bic_weightsd80),
                        pbias=rep(NA,8))
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfaicd90,dfaicd80,dfbic,dfbicd90,dfbicd80)
  
  return(dff)
  
}



stan_func_1 <- function(a,u) {
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                   S=dat$obsSpawners,
                   R=dat$obsRecruits,
                   logRS=log(dat$obsRecruits/dat$obsSpawners))
   #
  p <- samEst::ricker_stan(data=df,iter = 600, mod=m1f)
  #ricker autocorr
  pac <- samEst::ricker_stan(data=df,iter = 600, AC=TRUE, mod=m2f)
  #ricker tva
  ptva <- samEst::ricker_rw_stan(data=df, par="a",iter = 600, mod=m3f)
  #ricker tvab
  ptvb<- samEst::ricker_rw_stan(data=df, par="b",iter = 600, mod=m4f)
  
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",4)),each=nrow(df)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb"),each=nrow(df)),
                   by=rep(dat$year,4),
                   sim=rep(dat$alpha,4),
                   est=c(rep(p$alpha,nrow(df)),
                         rep(pac$alpha,nrow(df)),
                         ptva$alpha[2:c(nrow(df)+1)],
                         rep(mean(ptvb$alpha),nrow(df))))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",4)),each=nrow(df)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb"),each=nrow(df)),
                      by=rep(dat$year,4),
                      sim=rep(1/dat$beta,4),
                      est=c(rep(p$Smax,nrow(df)),
                            rep(pac$Smax,nrow(df)),
                            rep(ptva$Smax,nrow(df)),
                            ptvb$Smax))
  
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
                            
  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep(c(rep("MCMC",4)),each=nrow(df)),
                     model=rep(c("simple",
                                 "autocorr",
                                 "rwa","rwb"),each=nrow(df)),
                     by=rep(dat$year,4),
                     sim=rep(dat$sigma,4),
                     est=c(rep(p$sig,nrow(df)),
                           rep(pac$sig,nrow(df)),
                           rep(ptva$sig,nrow(df)),
                           rep(ptvb$sig,nrow(df))))
  
  dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
  #Pointwise loglikelihoods
  dfelpd<- data.frame(parameter="ELPD",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",4),
                     model=rep(c("simple",
                             "autocorr",
                             "rwa","rwb"),each=nrow(df)),
                     by=rep(dat$year,4),
                     sim=rep(NA,4),
                     est=c(
                       apply(rstan::extract(p$stanfit)$log_lik,2,samEst::log_mean_exp),
                       apply(rstan::extract(pac$stanfit)$log_lik,2,samEst::log_mean_exp),
                       apply(rstan::extract(ptva$stanfit)$log_lik,2,samEst::log_mean_exp),
                       apply(rstan::extract(ptvb$stanfit)$log_lik,2,samEst::log_mean_exp)),
                     pbias=rep(NA,4))
 
  
  dff<-rbind(dfa,dfsmax,dfsig,dfelpd)
  
  return(dff)
  
}

stan_func_2 <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                   S=dat$obsSpawners,
                   R=dat$obsRecruits,
                   logRS=log(dat$obsRecruits/dat$obsSpawners))
  
  #ricker tvab
  ptvab <- samEst::ricker_rw_stan(data=df, par="both",iter = 600, mod=m5f) 
  #ricker tvhmma
  phmma <- samEst::ricker_hmm_stan(data=df, par="a",iter = 600, mod=m6f)
  #ricker tvhmmb
  phmmb <- samEst::ricker_hmm_stan(data=df, par="b",iter = 600, mod=m7f)
  #ricker tvhmmab
  phmmab <- samEst::ricker_hmm_stan(data=df, par="both",iter = 600, mod=m8f) 
  
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",4)),each=nrow(df)),
                   model=rep(c("rwab",
                               "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                   by=rep(dat$year,4),
                   sim=rep(dat$alpha,4),
                   est=c(ptvab$alpha[2:c(nrow(df)+1)],
                         phmma$alpha[phmma$regime],
                         rep(phmmb$alpha,nrow(df)),
                         phmmab$alpha[phmm$regime]))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  dfsmax<- data.frame(parameter="Smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",4)),each=nrow(df)),
                      model=rep(c("rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                      by=rep(dat$year,4),
                      sim=rep(1/dat$beta,4),
                      est=c(ptvab$Smax,
                             rep(phmma$Smax,nrow(df)),
                             phmmb$Smax[phmmb$regime],
                             phmmab$Smax[phmm$regime]))
  
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep(c(rep("MCMC",4)),each=nrow(df)),
                     model=rep(c("rwab",
                                 "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                     by=rep(dat$year,4),
                     sim=rep(dat$sigma,4),
                     est=c(rep(ptvab$sigobs,nrow(df)),
                           rep(phmma$sigobs,nrow(df)),
                           rep(phmmb$sigobs,nrow(df)),
                           rep(phmmab$sigobs,nrow(df))))
  
  dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
  #Smsy - estimate from posterior 
  smsysim<-samEst::smsyCalc(dat$alpha,dat$beta)
  smsy_est=list()
  
  smsy_med1=NA
  for(t in 1:nrow(df)){smsy_med1[t]=median(samEst::smsyCalc(a=rstan::extract(ptvab$stanfit)$log_a[,t],b=rstan::extract(ptvab$stanfit)$b[,t]))}
  smsy_est[[1]]=smsy_med1
  #time-varying
  smsy_med2=NA
  smsy_med2[1]=median(samEst::smsyCalc(a=rstan::extract(phmma$stanfit)$log_a[,1],b=rstan::extract(phmma$stanfit)$b))
  smsy_med2[2]=median(samEst::smsyCalc(a=rstan::extract(phmma$stanfit)$log_a[,2],b=rstan::extract(phmma$stanfit)$b))
  
  smsy_est[[2]]=smsy_med2[phmma$regime]
  
  smsy_med3=NA
  smsy_med3[1]=median(samEst::smsyCalc(a=rstan::extract(phmmb$stanfit)$log_a,b=rstan::extract(phmmb$stanfit)$b[,1]))
  smsy_med3[2]=median(samEst::smsyCalc(a=rstan::extract(phmmb$stanfit)$log_a,b=rstan::extract(phmmb$stanfit)$b[,2]))
  
  smsy_est[[3]]=smsy_med3[phmmb$regime]
  
  smsy_med4=NA
  smsy_med4[1]=median(samEst::smsyCalc(a=rstan::extract(phmmab$stanfit)$log_a[,1],b=rstan::extract(phmmab$stanfit)$b[,1]))
  smsy_med4[2]=median(samEst::smsyCalc(a=rstan::extract(phmmab$stanfit)$log_a[,2],b=rstan::extract(phmmab$stanfit)$b[,2]))
  
  smsy_est[[4]]=smsy_med4[phmmab$regime]
  
  
  dfsmsy<- data.frame(parameter="smsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",4)),each=nrow(df)),
                      model=rep(c("rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                      by=rep(dat$year,4),
                      sim=rep(smsysim,4),
                      est=c(smsy_est[[1]],
                            smsy_est[[2]],
                            smsy_est[[3]],
                            smsy_est[[4]]))
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100
  
  
  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MCMC",4)),each=nrow(df)),
                       model=rep(c("rwab",
                                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                       by=rep(dat$year,4),
                       sim=rep(unlist(mapply(samEst::sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),4),
                       est=c(unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwab"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwab"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwab"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmma_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmma_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmb_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmab_regime"]))))
  
  dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100     
  #umsy
  
  ptvab$umsy=samEst::umsyCalc(unname(ptvab$alpha[-1]))
  phmma$umsy=samEst::umsyCalc(unname(phmma$alpha))
  phmmb$umsy=samEst::umsyCalc(unname(phmmb$alpha))
  phmmab$umsy=samEst::umsyCalc(unname(phmmab$alpha))
  
  dfumsy<- data.frame(parameter="umsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",4)),each=nrow(df)),
                      model=rep(c("rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                      by=rep(dat$year,4),
                      sim=rep(samEst::umsyCalc(dat$alpha),4),
                      est=c(ptvab$umsy,
                            phmma$umsy[phmma$regime],
                            rep(phmmb$umsy,nrow(df)),
                            phmmab$umsy[phmmab$regime])
                      )
  
  dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100
  
  #Pointwise loglikelihoods
  dfelpd<- data.frame(parameter="ELPD",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",4),
                      model=rep(c("rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
                      by=rep(dat$year,4),
                      sim=rep(NA,4),
                      est=c(
                        apply(rstan::extract(ptvab$stanfit)$log_lik,2,samEst::log_mean_exp),
                        apply(rstan::extract(phmma$stanfit)$log_lik,2,samEst::log_mean_exp),
                        apply(rstan::extract(phmmb$stanfit)$log_lik,2,samEst::log_mean_exp),
                        apply(rstan::extract(phmmab$stanfit)$log_lik,2,samEst::log_mean_exp)),
                      pbias=rep(NA,4))
  
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfelpd)
  
  return(dff)
  
}

pars<-data.frame(a=rep(1:3,each=2),
                 u=1:2)


m1f=samEst::sr_mod2(type='static',ac = FALSE,par='n',lfo =F)
m2f=samEst::sr_mod2(type='static',ac = TRUE,par='n',lfo=F)
m3f=samEst::sr_mod2(type='rw',par='a',lfo=F)
m4f=samEst::sr_mod2(type='rw',par='b',lfo=F)
m5f=samEst::sr_mod2(type='rw',par='both',lfo=F)
m6f=samEst::sr_mod2(type='hmm',par='a',lfo=F)
m7f=samEst::sr_mod2(type='hmm',par='b',lfo=F)
m8f=samEst::sr_mod2(type='hmm',par='both',lfo=F)



sjobstan1 <- slurm_apply(stan_func_1, pars, jobname = 'Stanrun',
                       nodes = 1, cpus_per_node = 12, submit = FALSE,
                       pkgs=c("samEst"),
                       rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                       libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                       global_objects=c("simPars","m1f","m2f","m3f","m4f"))

sjobstan2 <- slurm_apply(stan_func_1, pars, jobname = 'Stanrun2',
                        nodes = 1, cpus_per_node = 12, submit = FALSE,
                        pkgs=c("samEst"),
                        rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                        libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                        global_objects=c("simPars","m5f","m6f","m7f","m8f"))


save.image("./slrmjb.RData")

q()

#load again
library(rslurm)
setwd('..')
load('slrmjb.RData')

res3 <- get_slurm_out(sjobstan3, outtype = 'table', wait = FALSE)
head(res3, 3)


##cmdstan try

cmdstanr::set_cmdstan_path(path="/fs/vnas_Hdfo/comda/dag004/.cmdstan/cmdstan-2.31.0")
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1f.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f.stan")
mod3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5f.stan")
mod5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6f.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7f.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8f.stan")
mod8=cmdstanr::cmdstan_model(file8)


stan_func<- function(a,u){
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
            alpha_dirichlet=c(1,1)
            )
  
  #Clear temp. files from cmdstanr
  fl=list.files('/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  file.remove(paste('/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst/',fl,sep=''))
  
  #
  f1 <- mod1$sample(data=df,
                   seed = 123,
                   chains = 6, 
                   parallel_chains = 6,
                   iter_warmup = 200,
                   iter_sampling = 600,
                   refresh = 0,
                   adapt_delta = 0.95,
                   max_treedepth = 15)
  
  f2 <- mod2$sample(data=df,
                   seed = 123,
                   chains = 6, 
                   parallel_chains = 6,
                   iter_warmup = 200,
                   iter_sampling = 600,
                   refresh = 0,
                   adapt_delta = 0.95,
                   max_treedepth = 15)
  
  f3 <- mod3$sample(data=df,
                     seed = 123,
                     chains = 6, 
                     parallel_chains = 6,
                     iter_warmup = 200,
                     iter_sampling = 600,
                     refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f4 <- mod4$sample(data=df,
                      seed = 123,
                      chains = 6, 
                      parallel_chains = 6,
                      iter_warmup = 200,
                      iter_sampling = 600,
                      refresh = 0,
                      adapt_delta = 0.95,
                      max_treedepth = 15)
  
  f5 <- mod5$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f6 <- mod6$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 1,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f7 <- mod7$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  
  f8 <- mod8$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  #Max. prod
  
  #regime state sequence:
  phmma_alpha=f6$summary(variables=c('log_a'),'median')$median
  phmma_alpha_regime=phmma_alpha[f6$summary(variables=c('zstar'),'median')$median]
  phmmab_alpha=f8$summary(variables=c('log_a'),'median')$median
  phmmab_alpha_regime=phmmab_alpha[f8$summary(variables=c('zstar'),'median')$median]
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   est=c(rep(f1$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         rep(f2$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f3$summary(variables=c('log_a'),'median')$median,
                         rep(f4$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f5$summary(variables=c('log_a'),'median')$median,
                         phmma_alpha_regime,
                         rep(f7$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         phmmab_alpha_regime
                         ))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  #Smax
  phmmb_smax=f7$summary(variables=c('S_max'),'median')$median
  phmmb_smax_regime=phmmb_smax[f7$summary(variables=c('zstar'),'median')$median]
  phmmab_smax=f8$summary(variables=c('S_max'),'median')$median
  phmmab_smax_regime=phmmb_smax[f8$summary(variables=c('zstar'),'median')$median]
  
  dfsmax<- data.frame(parameter="smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(1/dat$beta,8),
                      est=c(rep(f1$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f2$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f3$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            f4$summary(variables=c('S_max'),'median')$median,
                            f5$summary(variables=c('S_max'),'median')$median,
                            rep(f6$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            phmmb_smax_regime,
                            phmmab_smax_regime
                            ))
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  #obs error
  dfsig<- data.frame(parameter="sigma_obs",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab",
                                 "hmma_regime","hmmb_regime","hmmab_regime"),
                     by=NA,
                     sim=NA,
                     est=c(f1$summary(variables=c('sigma'),'median')$median,
                           f2$summary(variables=c('sigma'),'median')$median,
                           f3$summary(variables=c('sigma'),'median')$median,
                           f4$summary(variables=c('sigma'),'median')$median,
                           f5$summary(variables=c('sigma'),'median')$median,
                           f6$summary(variables=c('sigma'),'median')$median,
                           f7$summary(variables=c('sigma'),'median')$median,
                           f8$summary(variables=c('sigma'),'median')$median))
  
  dfsig$pbias<- NA
  
  #sigma a
  dfsiga<- data.frame(parameter="sigma_a",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",2),
                     model=c("rwa","rwab"),
                     by=NA,
                     sim=NA,
                     est=c(f3$summary(variables=c('sigma_a'),'median')$median,
                           f5$summary(variables=c('sigma_a'),'median')$median))
  
  dfsiga$pbias<- NA
  
  #sigma b
  dfsigb<- data.frame(parameter="sigma_b",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwb","rwab"),
                      by=NA,
                      sim=NA,
                      est=c(f4$summary(variables=c('sigma_b'),'median')$median,
                            f5$summary(variables=c('sigma_b'),'median')$median))
  
  dfsigb$pbias<- NA
  #S msy
  #Smsy - estimate from posterior 
  #static
  smsysim<-samEst::smsyCalc(dat$alpha,dat$beta)

  dfsmsy<- data.frame(parameter="smsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(smsysim,8),
                      est=c(samEst::smsyCalc(a=dfa$est[dfa$model=="simple"],b=1/dfsmax$est[dfsmax$model=="simple"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="autocorr"],b=1/dfsmax$est[dfsmax$model=="autocorr"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwa"],b=1/dfsmax$est[dfsmax$model=="rwa"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwb"],b=1/dfsmax$est[dfsmax$model=="rwb"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwab"],b=1/dfsmax$est[dfsmax$model=="rwab"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmma_regime"],b=1/dfsmax$est[dfsmax$model=="hmma_regime"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmb_regime"],b=1/dfsmax$est[dfsmax$model=="hmmb_regime"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmab_regime"],b=1/dfsmax$est[dfsmax$model=="hmmab_regime"])))
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100
  
  
  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                       model=rep(c("simple",
                                   "autocorr",
                                   "rwa","rwb","rwab",
                                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                       by=rep(dat$year,8),
                       sim=rep(unlist(mapply(samEst::sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),8),
                       est=c(unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="simple"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="simple"], 
                                           b=1/dfsmax$est[dfsmax$model=="simple"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="autocorr"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="autocorr"], 
                                           b=1/dfsmax$est[dfsmax$model=="autocorr"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwa"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwa"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwa"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwb"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwb"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwb"])),
                         unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="rwab"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="rwab"], 
                                           b=1/dfsmax$est[dfsmax$model=="rwab"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmma_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmma_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmb_regime"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmab_regime"]))))
  
  dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100     
  
  #umsy
  dfumsy<- data.frame(parameter="umsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(samEst::umsyCalc(dat$alpha),8),
                      est=c(samEst::umsyCalc(dfa$est[dfa$model=="simple"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="autocorr"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwa"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwb"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwab"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmma_regime"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmb_regime"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmab_regime"]))
  )
  
  dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100
  
  #Pointwise loglikelihoods
  dfelpd<- data.frame(parameter="ELPD",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",8),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(NA,8),
                      est=c(
                        f1$summary(variables=c('log_lik'),'median')$median,
                        f2$summary(variables=c('log_lik'),'median')$median,
                        f3$summary(variables=c('log_lik'),'median')$median,
                        f4$summary(variables=c('log_lik'),'median')$median,
                        f5$summary(variables=c('log_lik'),'median')$median,
                        f6$summary(variables=c('log_lik'),'median')$median,
                        f7$summary(variables=c('log_lik'),'median')$median,
                        f8$summary(variables=c('log_lik'),'median')$median)
                        ,
                      pbias=rep(NA,8))
  
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsiga,dfsigb,dfsmsy,dfsgen,dfumsy,dfelpd)
  
  return(dff)
  
}


pars<-data.frame(a=rep(1:3,each=2),
                 u=1:2)


sjobstan3 <- slurm_apply(stan_func_3, pars, jobname = 'Stanrun3',
                         nodes = 1, cpus_per_node = 12, submit = FALSE,
                         pkgs=c("samEst","cmdstanr"),
                         rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                         libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                         global_objects=c("simPars","mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"))

save.image("./slrmjb.RData")

q()
#https://portal.science.gc.ca/confluence/display/SCIDOCS/Quick+Start+to+Using+Linux+Clusters+With+SLURM


#tester function
stan_func_t <- function(a,u) {
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
             alpha_dirichlet=c(1,1)
  )
  
  #Clear temp. files from cmdstanr
  unlink("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst/*")
  
  #
  f1 <- mod1$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f2 <- mod2$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f3 <- mod3$sample(data=df,
                    seed = 123,
                    
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f4 <- mod4$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f5 <- mod5$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f6 <- mod6$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 1,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f7 <- mod7$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  f8 <- mod8$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    parallel_chains = 6,
                    iter_warmup = 200,
                    iter_sampling = 600,
                    refresh = 1000,
                    max_treedepth = 20,
                    output_dir='/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst')
  
  p<- as.data.frame(f1$draws(format='df'))
  pac<- as.data.frame(f2$draws(format='df'))
  ptva<- as.data.frame(f3$draws(format='df'))
  ptvb<- as.data.frame(f4$draws(format='df'))
  ptvab<- as.data.frame(f5$draws(format='df'))
  phmma<- as.data.frame(f6$draws(format='df'))
  phmmb<- as.data.frame(f7$draws(format='df'))
  phmmab<- as.data.frame(f8$draws(format='df'))
  
  #Max. prod
  ptva_alpha=apply(ptva[grepl('log_a',colnames(ptva))],2,median)[-1]
  ptvab_alpha=apply(ptvab[grepl('log_a',colnames(ptvab))],2,median)[-1]
  phmma_alpha=apply(phmma[grepl('log_a',colnames(phmma))],2,median)
  phmma_alpha_regime=phmma_alpha[apply(phmma[grepl('zstar',colnames(phmma))],2,median)[-c(nrow(dat)+1)]]
  phmmab_alpha=apply(phmmab[grepl('log_a',colnames(phmmab))],2,median)
  phmmab_alpha_regime=phmmab_alpha[apply(phmmab[grepl('zstar',colnames(phmmab))],2,median)[-c(nrow(dat)+1)]]
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(dat)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   est=c(rep(median(p$log_a),nrow(dat)),
                         rep(median(pac$log_a),nrow(dat)),
                         ptva_alpha,
                         rep(median(ptvb$log_a),nrow(dat)),
                         ptvab_alpha,
                         phmma_alpha_regime,
                         rep(median(phmmb$log_a),nrow(dat)),
                         phmmab_alpha_regime
                   ))
  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
  
  #Clear temp. files from cmdstanr
  unlink("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp cmdst/*")
  
  return(dfa)
  
}


sjobstan3 <- slurm_apply(stan_func_t, pars, jobname = 'Stanrun3',
                         nodes = 1, cpus_per_node = 12, submit = FALSE,
                         pkgs=c("samEst"),
                         rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                         libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                         global_objects=c("simPars","mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8"))



###LOCAL pbias & AIC/BIC"####
#test it out
tmb_func <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfaic,dfbic)
  
  return(dff)
  
}

pars<-data.frame(a=rep(1:12,each=5000),
                 u=1:5000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                       nodes = 1, cpus_per_node = 50, submit = FALSE,
                       pkgs=c("samEst"),
                       rscript_path = "/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/",
                       libPaths="/fs/vnas_Hdfo/comda/dag004/Rlib/",
                       global_objects=c("simPars"))

save.image("./slrmjb.RData")


####TMB LFO FUNC####
tmb_func <- function(a,u) {
  
  allsimest <- list()
  simData<- readRDS(paste0("/fs/vnas_Hdfo/comda/dag004/Documents/cluster-tvsimest/outs/SamSimOutputs/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
  lfoalpha <- tryCatch(tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=10),error = function(e) {lfoalpha=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last5param=rep(-999,length(lfoac$lastparam)))})
  lfobeta <- tryCatch(tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=10),error = function(e) {lfobeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last5param=rep(-999,length(lfoac$lastparam)))})
  lfoalphabeta <- tryCatch(tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=10),error = function(e) {lfoalphabeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last5param=rep(-999,length(lfoac$lastparam)))})
  lfohmma <- tryCatch(tmb_mod_lfo_cv(data=df,model='HMM_a', L=10),error = function(e) {lfohmma=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last5regime_pick=rep(-999,length(lfoac$lastparam)))})
  lfohmmb <- tryCatch(tmb_mod_lfo_cv(data=df,model='HMM_b', L=10),error = function(e) {lfohmmb=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last5regime_pick=rep(-999,length(lfoac$lastparam)))})
  lfohmm <- tryCatch(tmb_mod_lfo_cv(data=df,model='HMM', L=10),error = function(e) {lfohmm=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                last5regime_pick=rep(-999,length(lfoac$lastparam)))})
}
  p <- ricker_TMB(data=df)
  pac <- ricker_TMB(data=df, AC=TRUE)
  ptva <- ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- ricker_rw_TMB(data=df, tv.par='both')
  phmma <- ricker_hmm_TMB(data=df, tv.par='a')
  phmmb <- ricker_hmm_TMB(data=df, tv.par='b')
  phmm  <- ricker_hmm_TMB(data=df, tv.par='both')
  
  
  
  
  
  dso_filename = m2f@dso@dso_filename
  loaded_dlls = getLoadedDLLs()
  if (dso_filename %in% names(loaded_dlls)) {
    message("Unloading DLL for model dso ", dso_filename)
    model.dll = loaded_dlls[[dso_filename]][['path']]
    dyn.unloam1fd(model.dll)
  } else {
    message("No loaded DLL for model dso ", dso_filename)
  }
  
  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[stringr::str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  message("DLL Count = ", length(getLoadedDLLs()), ": [", stringr::str_c(names(loaded_dlls), collapse = ","), "]")
  