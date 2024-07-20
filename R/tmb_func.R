
tmb_func <- function(path=".",a, u) {
  
  print(paste("a is ",a,"u is ",u))
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

  Smax_mean<-(max(df$S)*.5)
  Smax_sd<-Smax_mean
 
  logbeta_pr_sig = sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_pr = log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
  
  dirpr <- matrix(c(2,1,1,2),2,2)


  p <- tryCatch({ ricker_TMB(data=df,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  pac <- tryCatch({ricker_TMB(data=df, AC=TRUE,logb_p_mean=logbeta_pr,
                 logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptva <- tryCatch({ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,
                  logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvb <- tryCatch({ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})
  
  ptvab <- tryCatch({ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmma <- tryCatch({ricker_hmm_TMB(data=df, tv.par='a', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmmb <- tryCatch({ricker_hmm_TMB(data=df, tv.par='b', dirichlet_prior=dirpr,
                    logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmm <- tryCatch({ricker_hmm_TMB(data=df, tv.par='both', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )
  
  dfa <- data.frame(parameter="logalpha",
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
              mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$logalpha}, nrow(df)),
                    rep(if(!is.null(pac$fail_conv)){NA}else{pac$logalpha}, nrow(df)),
                    if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$logalpha},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$logalpha},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$logalpha},
                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$logalpha[phmma$regime]},
                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$logalpha}, nrow(df)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$logalpha[phmm$regime]}), 
              convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
              conv_warning=rep(c( p$conv_problem,
                    pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                    phmma$conv_problem,
                    phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax <- data.frame(parameter="Smax",
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
      mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$Smax}, nrow(df)),
        rep(if(!is.null(pac$fail_conv)){NA}else{pac$Smax}, nrow(df)),
        if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smax},
        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smax},
        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smax},
        rep(if(!is.null(phmma$fail_conv)){NA}else{phmma$Smax}, nrow(df)),
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smax[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smax[phmm$regime]}),
      convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
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
      mode=rep(c(ifelse(is.null(p$fail_conv),p$sig,NA),
                 ifelse(is.null(pac$fail_conv),pac$sig,NA),
                 ifelse(is.null(ptva$fail_conv),ptva$sig,NA),
                 ifelse(is.null(ptvb$fail_conv),ptvb$sig,NA),
                 ifelse(is.null(ptvab$fail_conv),ptvab$sig,NA),
                 ifelse(is.null(phmma$fail_conv),phmma$sigma,NA),
                 ifelse(is.null(phmmb$fail_conv),phmmb$sigma,NA),
                 ifelse(is.null(phmm$fail_conv),phmm$sigma,NA)),each=nrow(df)), 
      convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df)))
    
    dfsig$pbias <- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias <- (dfsig$mode-dfsig$sim)

    #sigma_a
    dfsiga<- data.frame(parameter="sigma_a",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=c("rwa","rwab"),
      by=NA,
      sim=NA,
      median=NA,
      mode=c(ifelse(is.null(ptva$fail_conv),ptva$siga,NA),
        ifelse(is.null(ptvab$fail_conv),ptvab$siga,NA)),
      convergence=c( ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv)),
      conv_warning=c( ptva$conv_problem,
         ptvab$conv_problem))
    
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
      mode=c(ifelse(is.null(ptvb$fail_conv),ptvb$sigb,NA),
        ifelse(is.null(ptvab$fail_conv),ptvab$sigb,NA)),
      convergence=c( ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv)),
      conv_warning=c(ptvb$conv_problem,        
        ptvab$conv_problem)
      )
    
    dfsigb$pbias <- ((dfsigb$mode-dfsigb$sim)/dfsigb$sim)*100
    dfsigb$bias <- (dfsigb$mode-dfsigb$sim)
              
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
      mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$Smsy}, nrow(df)),
        rep(if(!is.null(pac$fail_conv)){NA}else{pac$Smsy}, nrow(df)),
        if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smsy},
        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smsy},
        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smsy},
        if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$Smsy[phmma$regime]},
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smsy[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smsy[phmm$regime]}),    
        convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
        conv_warning=rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df))) 
  
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
    mode=c(
      if(is.null(p$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="simple"&dfa$method=="MLE"],
            Smsy=dfsmsy$mode[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
            b=1/dfsmax$mode[dfsmax$model=="simple"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},
      
      if(is.null(pac$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="autocorr"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="autocorr"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(ptva$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwa"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(ptvb$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(ptvab$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmma$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmmb$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmm$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))}),
     
    convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
    conv_warning=rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df)))
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
 
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
    mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$umsy}, nrow(df)),
                    rep(if(!is.null(pac$fail_conv)){NA}else{pac$umsy}, nrow(df)),
                    if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$umsy},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$umsy},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$umsy},
                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$umsy[phmma$regime]},
                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$umsy}, nrow(df)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$umsy[phmm$regime]}), 
    convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv)
                    ),each=nrow(df)),
    conv_warning=rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df)))

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
                       mode=c(ifelse(is.null(p$fail_conv),p$AICc, NA),
                              ifelse(is.null(pac$fail_conv), pac$AICc, NA),
                              ifelse(is.null(ptva$fail_conv), ptva$AICc, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvb$AICc, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvab$AICc, NA),
                              ifelse(is.null(phmma$fail_conv), phmma$AICc, NA),
                              ifelse(is.null(phmmb$fail_conv), phmmb$AICc, NA),
                              ifelse(is.null(phmm$fail_conv), phmm$AICc, NA)),
                       convergence=c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                              ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                              ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                              ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                              ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                              ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                              ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                              ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv) ),
                       conv_warning= c( p$conv_problem,
                                      pac$conv_problem,
                                      ptva$conv_problem,
                                      ptvb$conv_problem,
                                      ptvab$conv_problem,
                                      phmma$conv_problem,
                                      phmmb$conv_problem,
                                      phmm$conv_problem),
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
                       mode=c(ifelse(is.null(p$fail_conv),p$BIC, NA),
                              ifelse(is.null(pac$fail_conv), pac$BIC, NA),
                              ifelse(is.null(ptva$fail_conv), ptva$BIC, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvb$BIC, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvab$BIC, NA),
                              ifelse(is.null(phmma$fail_conv), phmma$BIC, NA),
                              ifelse(is.null(phmmb$fail_conv), phmmb$BIC, NA),
                              ifelse(is.null(phmm$fail_conv), phmm$BIC, NA)),
                       convergence=c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                                     ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                                     ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                                     ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                                     ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                                     ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                                     ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                                     ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv) ),
                       conv_warning= c( p$conv_problem,
                                      pac$conv_problem,
                                      ptva$conv_problem,
                                      ptvb$conv_problem,
                                      ptvab$conv_problem,
                                      phmma$conv_problem,
                                      phmmb$conv_problem,
                                      phmm$conv_problem),
                       pbias=rep(NA,8),
                       bias=rep(NA,8))

    #EDF
    dfedf<- data.frame(parameter="EDF",
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
                       mode=c(ifelse(is.null(p$fail_conv),3, NA),
                              ifelse(is.null(pac$fail_conv),4, NA),
                              ifelse(is.null(ptva$fail_conv), ptva$EDF, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvb$EDF, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvab$EDF, NA),
                              ifelse(is.null(phmma$fail_conv), length(unlist(phmma$tmb_params)), NA),
                              ifelse(is.null(phmmb$fail_conv), length(unlist(phmmb$tmb_params)), NA),
                              ifelse(is.null(phmm$fail_conv), length(unlist(phmm$tmb_params)), NA)),
                       convergence=c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                                     ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                                     ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                                     ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                                     ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                                     ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                                     ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                                     ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv) ),
                       conv_warning= c( p$conv_problem,
                                      pac$conv_problem,
                                      ptva$conv_problem,
                                      ptvb$conv_problem,
                                      ptvab$conv_problem,
                                      phmma$conv_problem,
                                      phmmb$conv_problem,
                                      phmm$conv_problem),
                       pbias=rep(NA,8),
                       bias=rep(NA,8))
  

  
   #lfo
    lfostatic <- tmb_mod_lfo_cv(data=df,model='static', L=15, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfoac <- tmb_mod_lfo_cv(data=df,model='staticAC', L=15, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfoalpha <- tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=15,priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfobeta <- tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=15, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfoalphabeta <- tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=15, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfohmma <- tmb_mod_lfo_cv(data=df,model='HMM_a', L=15, dirichlet_prior=dirpr, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfohmmb <- tmb_mod_lfo_cv(data=df,model='HMM_b', L=15, dirichlet_prior=dirpr, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    lfohmm <- tmb_mod_lfo_cv(data=df,model='HMM', L=15, dirichlet_prior=dirpr, priorlogb="default",logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
    
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
                           sum(lfohmm$conv_problem) ),
                       conv_warning=NA,
                       pbias=rep(NA,20),
                       bias=rep(NA,20))

    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfsiga,dfsigb,
      dfaic,dfbic,dflfo,dfedf)

  return(dff)

}
