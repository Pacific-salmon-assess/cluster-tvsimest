
tmb_func_retro <- function(path=".",a, u, minyr=10) {
  
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

  
  dirpr <- matrix(c(2,1,1,2),2,2)


  dfaretro<-list()
  dfsmaxretro<-list()
  dfsigretro<-list()
  dfsigaretro<-list()
  dfsigbretro<-list()
  dfsmsyretro<-list()
  dfsgenretro<-list()
  dfumsyretro<-list()
  dfaicretro<-list()
  dfbicretro<-list()

  smsysim<-smsyCalc(dat$alpha,dat$beta)

  for(i in minyr:nrow(df)){
    dfset<-df[1:i,]
    ct<-i-minyr+1
    
    Smax_mean<-(max(dfset$S)*.5)
    Smax_sd<-Smax_mean
 
    logbeta_pr_sig = sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
    logbeta_pr = log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
  
    pac <- tryCatch({ricker_TMB(data=dfset, AC=TRUE,logb_p_mean=logbeta_pr,
                                logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                    error=function(cond){
                      message(cond)
                      return(list(fail_conv=1,
                                  conv_problem=1))})

    ptva <- tryCatch({ricker_rw_TMB(data=dfset,tv.par='a',logb_p_mean=logbeta_pr,
                  logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

    ptvb <- tryCatch({ricker_rw_TMB(data=dfset, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})
  
    ptvab <- tryCatch({ricker_rw_TMB(data=dfset, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

    phmma <- tryCatch({ricker_hmm_TMB(data=dfset, tv.par='a', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

    phmmb <- tryCatch({ricker_hmm_TMB(data=dfset, tv.par='b', dirichlet_prior=dirpr,
                    logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

    phmm <- tryCatch({ricker_hmm_TMB(data=dfset, tv.par='both', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )
    
    dfa <- data.frame(parameter="logalpha",
              endyr=i,
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",7)),each=nrow(dfset)),
              model=rep(c('autocorr',"rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=nrow(dfset)),
              by=rep(dat$year[1:i],7),
              sim=rep(dat$alpha[1:i],7),
              median=NA,
              mode=c(if(!is.null(pac$fail_conv)){rep(NA, nrow(dfset))}else{pac$logalpha}
                    if(!is.null(ptva$fail_conv)){rep(NA, nrow(dfset))}else{ptva$logalpha},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(dfset))}else{ptvb$logalpha},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(dfset))}else{ptvab$logalpha},
                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(dfset))}else{phmma$logalpha[phmma$regime]},
                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$logalpha}, nrow(dfset)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(dfset))}else{phmm$logalpha[phmm$regime]}), 
              convergence=rep(c(
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(dfset)),
              conv_warning=rep(c( 
                    pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                    phmma$conv_problem,
                    phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(dfset)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax <- data.frame(parameter="Smax",
      endyr=i,
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",7)),each=nrow(dfset)),
      model=rep(c('autocorr',"rwa","rwb","rwab",
        "hmma","hmmb","hmmab"),each=nrow(dfset)),
      by=rep(dat$year[1:i],7),
      sim=rep(1/dat$beta[1:i],7),
      median=NA,
      mode=c(
        if(!is.null(pac$fail_conv)){rep(NA, nrow(dfset))}else{pac$Smax},
        if(!is.null(ptva$fail_conv)){rep(NA, nrow(dfset))}else{ptva$Smax},
        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(dfset))}else{ptvb$Smax},
        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(dfset))}else{ptvab$Smax},
        rep(if(!is.null(phmma$fail_conv)){NA}else{phmma$Smax}, nrow(dfset)),
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(dfset))}else{phmmb$Smax[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(dfset))}else{phmm$Smax[phmm$regime]}),
      convergence=rep(c(
        ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
        ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(dfset)),
      conv_warning=rep(c(
        pac$conv_problem,
        ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(dfset)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
   
    #Smsy
    
  
    dfsmsy<- data.frame(parameter="smsy",
      endyr=i,
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",7)),each=nrow(dfset)),
      model=rep(c('autocorr', "rwa","rwb","rwab",
        "hmma","hmmb","hmmab"),each=nrow(dfset)),
      by=rep(dat$year[1:i],7),
      sim=rep(smsysim[1:i],7),
      median=NA,
      mode=c(
        if(!is.null(pac$fail_conv)){rep(NA, nrow(dfset))}else{pac$Smsy},
        if(!is.null(ptva$fail_conv)){rep(NA, nrow(dfset))}else{ptva$Smsy},
        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(dfset))}else{ptvb$Smsy},
        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(dfset))}else{ptvab$Smsy},
        if(!is.null(phmma$fail_conv)){rep(NA, nrow(dfset))}else{phmma$Smsy[phmma$regime]},
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(dfset))}else{phmmb$Smsy[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(dfset))}else{phmm$Smsy[phmm$regime]}),    
        convergence=rep(
                  c(ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(dfset)),
        conv_warning=rep(c(
                    pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(dfset))) 
  
  dfsmsy$pbias<- ((dfsmsy$mode-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (dfsmsy$mode-dfsmsy$sim)

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
    endyr=i,
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",7)),each=nrow(dfset)),
    model=rep(c('autocorr',"rwa","rwb","rwab",
      "hmma","hmmb","hmmab"),each=nrow(dfset)),
    by=rep(dat$year[1:i],7),
    sim=rep(unlist(mapply(sGenCalc,a=dat$alpha[1:i],Smsy=smsysim[1:i], b=dat$beta[1:i])),7),
    median=NA,
    mode=c(
      if(is.null(pac$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="autocorr"&dfa$method=="MLE"],
                                                Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
                                                b=1/dfsmax$mode[dfsmax$model=="autocorr"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},
      
      
          if(is.null(ptva$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwa"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},

       if(is.null(ptvb$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},

       if(is.null(ptvab$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},

       if(is.null(phmma$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},

       if(is.null(phmmb$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))},

       if(is.null(phmm$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(dfset))}),
     
    convergence=rep(c(
                   ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(dfset)),
    conv_warning=rep(c(pac$conv_problem, 
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(dfset)))
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
 
  dfumsy<- data.frame(parameter="umsy",
    endyr=i,
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",7)),each=nrow(dfset)),
    model=rep(c('autocorr',"rwa","rwb","rwab",
      "hmma","hmmb","hmmab"),each=nrow(dfset)),
    by=rep(dat$year[1:i],7),
    sim=rep(umsyCalc(dat$alpha)[1:i],7),
    median=NA,
    mode=c(
      if(!is.null(pac$fail_conv)){rep(NA, nrow(dfset))}else{pac$umsy},
      if(!is.null(ptva$fail_conv)){rep(NA, nrow(dfset))}else{ptva$umsy},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(dfset))}else{ptvb$umsy},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(dfset))}else{ptvab$umsy},
                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(dfset))}else{phmma$umsy[phmma$regime]},
                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$umsy}, nrow(dfset)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(dfset))}else{phmm$umsy[phmm$regime]}), 
    convergence=rep(c(
                 ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(dfset)),
    conv_warning=rep(c(
                    pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(dfset)))

    dfumsy$pbias<- ((dfumsy$mode-dfumsy$sim)/dfumsy$sim)*100
    dfumsy$bias<- (dfumsy$mode-dfumsy$sim)

   
    
    dfaretro[[ct]] <- dfa
    dfsmaxretro[[ct]] <- dfsmax
    dfsmsyretro[[ct]] <- dfsmsy
    dfsgenretro[[ct]] <- dfsgen
    dfumsyretro[[ct]] <- dfumsy
 

  }

  rtr_loga <- do.call(rbind,dfaretro)
  rtr_smax <- do.call(rbind,dfsmaxretro)
  rtr_smsy <- do.call(rbind,dfsmsyretro)
  rtr_sgen <- do.call(rbind,dfsgenretro)
  rtr_umsy <- do.call(rbind,dfumsyretro)


  dff<-rbind(rtr_loga,rtr_smax,rtr_smsy,
    rtr_sgen,rtr_umsy)

  return(dff)

}
