
tmb_func_comp_hmm_smax <- function(path=".",a, u) {
  
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


  phmma2 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='a', dirichlet_prior=dirpr,
                  Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmmb2 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='b', dirichlet_prior=dirpr,
                     Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmm2 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='both', dirichlet_prior=dirpr,
                   Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )
  

  phmma3 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='a', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmmb3 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='b', dirichlet_prior=dirpr,
                    logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmm3 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='both', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )

  dfa <- data.frame(parameter="logalpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",9)),each=nrow(df)),
              model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
              by=rep(dat$year,9),
              sim=rep(dat$alpha,9),
              median=NA,
              mode=c(if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmma$logalpha[phmma$regime]},
                   rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$logalpha}, nrow(df)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$logalpha[phmm$regime]},
                   if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$logalpha[phmma2$regime]},
                    rep(if(!is.null(phmmb2$fail_conv)){NA}else{phmmb2$logalpha}, nrow(df)),
                    if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$logalpha[phmm2$regime]},
                   if(!is.null(phmma3$fail_conv)){rep(NA, nrow(df))}else{phmma3$logalpha[phmma3$regime]},
                    rep(if(!is.null(phmmb3$fail_conv)){NA}else{phmmb3$logalpha}, nrow(df)),
                    if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$logalpha[phmm3$regime]}
                  ), 
              convergence=rep(c(
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
              conv_warning=rep(c(phmma$conv_problem,
                    phmmb$conv_problem,
                    phmm$conv_problem,
                    phmma2$conv_problem,
                    phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                    phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax <- data.frame(parameter="Smax",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",9)),each=nrow(df)),
      model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
      by=rep(dat$year,9),
      sim=rep(1/dat$beta,9),
      median=NA,
      mode=c(
        rep(if(!is.null(phmma$fail_conv)){NA}else{phmma$Smax}, nrow(df)),
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smax[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smax[phmm$regime]},
        rep(if(!is.null(phmma2$fail_conv)){NA}else{phmma2$Smax}, nrow(df)),
        if(!is.null(phmmb2$fail_conv)){rep(NA, nrow(df))}else{phmmb2$Smax[phmmb2$regime]},
        if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$Smax[phmm2$regime]},
        rep(if(!is.null(phmma3$fail_conv)){NA}else{phmma3$Smax}, nrow(df)),
        if(!is.null(phmmb3$fail_conv)){rep(NA, nrow(df))}else{phmmb3$Smax[phmmb3$regime]},
        if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$Smax[phmm3$regime]}
      ),
      convergence=rep(c(ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c(phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem,
                    phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
      by=rep(dat$year,9),
      sim=rep(dat$sigma,9),
      median=NA,
      mode=rep(c(ifelse(is.null(phmma$fail_conv),phmma$sigma,NA),
                 ifelse(is.null(phmmb$fail_conv),phmmb$sigma,NA),
                 ifelse(is.null(phmm$fail_conv),phmm$sigma,NA),
                 ifelse(is.null(phmma2$fail_conv),phmma2$sigma,NA),
                 ifelse(is.null(phmmb2$fail_conv),phmmb2$sigma,NA),
                 ifelse(is.null(phmm2$fail_conv),phmm2$sigma,NA),
                 ifelse(is.null(phmma3$fail_conv),phmma2$sigma,NA),
                 ifelse(is.null(phmmb3$fail_conv),phmmb2$sigma,NA),
                 ifelse(is.null(phmm3$fail_conv),phmm2$sigma,NA)
               ),each=nrow(df)), 
      convergence=rep(c(ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c(phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem,
                     phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
    
    dfsig$pbias <- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias <- (dfsig$mode-dfsig$sim)

              
    #Smsy
    smsysim<-smsyCalc(dat$alpha,dat$beta)
  
    dfsmsy<- data.frame(parameter="smsy",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",9)),each=nrow(df)),
      model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
      by=rep(dat$year,9),
      sim=rep(smsysim,9),
      median=NA,
      mode=c(if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$Smsy[phmma$regime]},
        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smsy[phmmb$regime]},
        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smsy[phmm$regime]},
        if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$Smsy[phmma2$regime]},
        if(!is.null(phmmb2$fail_conv)){rep(NA, nrow(df))}else{phmmb2$Smsy[phmmb2$regime]},
        if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$Smsy[phmm2$regime]},
        if(!is.null(phmma3$fail_conv)){rep(NA, nrow(df))}else{phmma3$Smsy[phmma3$regime]},
        if(!is.null(phmmb3$fail_conv)){rep(NA, nrow(df))}else{phmmb3$Smsy[phmmb3$regime]},
        if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$Smsy[phmm3$regime]}
      ),    
        convergence=rep(c(ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
        conv_warning=rep(c(phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem,
                     phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df))) 
  
  dfsmsy$pbias<- ((dfsmsy$mode-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (dfsmsy$mode-dfsmsy$sim)

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",9)),each=nrow(df)),
    model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
    by=rep(dat$year,9),
    sim=rep(unlist(mapply(sGenCalc,loga=dat$alpha,Smsy=smsysim, b=dat$beta)),9),
    median=NA,
    mode=c(
      if(is.null(phmma$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmmb$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmm$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmma2$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmma_all_smax"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma_all_smax"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma_all_smax"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmmb2$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmb_all_smax"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb_all_smax"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb_all_smax"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmm2$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmab_all_smax"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab_all_smax"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab_all_smax"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},
       
       if(is.null(phmma3$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmma_all"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma_all"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmma_all"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmmb3$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmb_all"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb_all"&dfsmsy$method=="MLE"],
           b=1/dfsmax$mode[dfsmax$model=="hmmb_all"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))},

       if(is.null(phmm3$fail_conv)){unlist(mapply(sGenCalc,loga=dfa$mode[dfa$model=="hmmab_all"&dfa$method=="MLE"],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab_all"&dfsmsy$method=="MLE"], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab_all"&dfsmax$method=="MLE"]))}else{rep(NA, nrow(df))}
     ),
     
    convergence=rep(c(
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
    conv_warning=rep(c( phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem,
                    phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                     phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
 
  dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",9)),each=nrow(df)),
    model=rep(c( "hmma","hmmb","hmmab",
                   "hmma_all_smax","hmmb_all_smax","hmmab_all_smax",
                    "hmma_all","hmmb_all","hmmab_all"),each=nrow(df)),
    by=rep(dat$year,9),
    sim=rep(umsyCalc(dat$alpha),9),
    median=NA,
    mode=c(
                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$umsy[phmma$regime]},
                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$umsy}, nrow(df)),
                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$umsy[phmm$regime]},
                    if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$umsy[phmma2$regime]},
                    rep(if(!is.null(phmmb2$fail_conv)){NA}else{phmmb2$umsy}, nrow(df)),
                    if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$umsy[phmm2$regime]},
                    if(!is.null(phmma3$fail_conv)){rep(NA, nrow(df))}else{phmma3$umsy[phmma3$regime]},
                    rep(if(!is.null(phmmb3$fail_conv)){NA}else{phmmb3$umsy}, nrow(df)),
                    if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$umsy[phmm3$regime]}
                  ), 
    convergence=rep(c( ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
    conv_warning=rep(c( phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem,
                    phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))

    dfumsy$pbias<- ((dfumsy$mode-dfumsy$sim)/dfumsy$sim)*100
    dfumsy$bias<- (dfumsy$mode-dfumsy$sim)

    
    


    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy)

  return(dff)

}
