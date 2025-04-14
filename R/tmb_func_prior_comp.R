
tmb_func_prior_comp <- function(path=".",a, u) {
  
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

  p2 <- tryCatch({ ricker_TMB(data=df, priors_flag=0, silent=TRUE)},
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

  pac2 <- tryCatch({ricker_TMB(data=df, AC=TRUE, priors_flag=0, silent=TRUE)},
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

  ptva2 <- tryCatch({ricker_rw_TMB(data=df,tv.par='a', priors_flag=0, deltaEDF=0.0001, silent=TRUE)},
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

  ptvb2 <- tryCatch({ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
                    priors_flag=0, deltaEDF=0.0001, silent=TRUE)},
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

  ptvab2 <- tryCatch({ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,
                   priors_flag=0, deltaEDF=0.0001, silent=TRUE)},
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

  phmma2 <- tryCatch({ricker_hmm_TMB(data=df, tv.par='a', dirichlet_prior=dirpr,
                  priors_flag=0, silent=TRUE)},
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

  phmmb2 <- tryCatch({ricker_hmm_TMB(data=df, tv.par='b', dirichlet_prior=dirpr,
                    priors_flag=0, silent=TRUE)},
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

  phmm2 <- tryCatch({ricker_hmm_TMB(data=df, tv.par='both', dirichlet_prior=dirpr,
                  priors_flag=0, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )
  
  dfa <- data.frame(parameter="logalpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",16)),each=nrow(df)),
              model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
              prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
              by=rep(dat$year,16),
              sim=rep(dat$alpha,16),
              median=NA,
              mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$logalpha}, nrow(df)),
                    rep(if(!is.null(p2$fail_conv)){NA}else{p2$logalpha}, nrow(df)),
                    
                    rep(if(!is.null(pac$fail_conv)){NA}else{pac$logalpha}, nrow(df)),
                    rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$logalpha}, nrow(df)),
                    
                    if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$logalpha},
                    if(!is.null(ptva2$fail_conv)){rep(NA, nrow(df))}else{ptva2$logalpha},

                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$logalpha},
                    if(!is.null(ptvb2$fail_conv)){rep(NA, nrow(df))}else{ptvb2$logalpha},

                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$logalpha},
                    if(!is.null(ptvab2$fail_conv)){rep(NA, nrow(df))}else{ptvab2$logalpha},

                    if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$logalpha[phmma$regime]},
                    if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$logalpha[phmma$regime]},

                    rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$logalpha}, nrow(df)),
                    rep(if(!is.null(phmmb2$fail_conv)){NA}else{phmmb2$logalpha}, nrow(df)),

                    if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$logalpha[phmm$regime]},
                    if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$logalpha[phmm$regime]}

                  ), 
              convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv)

                    ),each=nrow(df)),
              conv_warning=rep(c( p$conv_problem,p2$conv_problem,
                    pac$conv_problem,  pac2$conv_problem,  
                    ptva$conv_problem, ptva2$conv_problem,  
                    ptvb$conv_problem, ptvb2$conv_problem,  
                    ptvab$conv_problem,ptvab2$conv_problem,
                    phmma$conv_problem,phmma2$conv_problem,
                    phmmb$conv_problem,phmmb2$conv_problem,
                    phmm$conv_problem,  phmm2$conv_problem
                    ),each=nrow(df)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax <- data.frame(parameter="Smax",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",16)),each=nrow(df)),
      model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
      prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
      by=rep(dat$year,16),
      sim=rep(1/dat$beta,16),
      median=NA,
      mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$Smax}, nrow(df)), 
             rep(if(!is.null(p2$fail_conv)){NA}else{p2$Smax}, nrow(df)),

        rep(if(!is.null(pac$fail_conv)){NA}else{pac$Smax}, nrow(df)), 
        rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$Smax}, nrow(df)),

        if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smax},
        if(!is.null(ptva2$fail_conv)){rep(NA, nrow(df))}else{ptva2$Smax},

        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smax},
        if(!is.null(ptvb2$fail_conv)){rep(NA, nrow(df))}else{ptvb2$Smax},

        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smax},
        if(!is.null(ptvab2$fail_conv)){rep(NA, nrow(df))}else{ptvab2$Smax},

        rep(if(!is.null(phmma$fail_conv)){NA}else{phmma$Smax}, nrow(df)),
        rep(if(!is.null(phmma2$fail_conv)){NA}else{phmma2$Smax}, nrow(df)),

        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smax[phmmb$regime]},
        if(!is.null(phmmb2$fail_conv)){rep(NA, nrow(df))}else{phmmb2$Smax[phmmb$regime]},

        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smax[phmm$regime]},
        if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$Smax[phmm$regime]}
      ),
      convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                       ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),
                    
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),
                    
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv), 
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv), 

                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c( p$conv_problem,p2$conv_problem,
                     pac$conv_problem,   pac2$conv_problem,   
                     ptva$conv_problem,  ptva2$conv_problem,   
                     ptvb$conv_problem,  ptvb2$conv_problem,   
                     ptvab$conv_problem, ptvab2$conv_problem,  
                     phmma$conv_problem, phmma2$conv_problem,   
                     phmmb$conv_problem, phmmb2$conv_problem,   
                    phmm$conv_problem , phmm2$conv_problem    
                    ),each=nrow(df)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
      prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
      by=rep(dat$year,16),
      sim=rep(dat$sigma,16),
      median=NA,
      mode=rep(c(ifelse(is.null(p$fail_conv),p$sig,NA),
                 ifelse(is.null(p2$fail_conv),p2$sig,NA),
                 ifelse(is.null(pac$fail_conv),pac$sig,NA),
                 ifelse(is.null(pac2$fail_conv),pac2$sig,NA),
                 ifelse(is.null(ptva$fail_conv),ptva$sig,NA),
                 ifelse(is.null(ptva2$fail_conv),ptva2$sig,NA),
                 ifelse(is.null(ptvb$fail_conv),ptvb$sig,NA),
                 ifelse(is.null(ptvb2$fail_conv),ptvb2$sig,NA),
                 ifelse(is.null(ptvab$fail_conv),ptvab$sig,NA),
                 ifelse(is.null(ptvab2$fail_conv),ptvab2$sig,NA),
                 ifelse(is.null(phmma$fail_conv),phmma$sigma,NA),
                 ifelse(is.null(phmma2$fail_conv),phmma2$sigma,NA),
                 ifelse(is.null(phmmb$fail_conv),phmmb$sigma,NA),
                 ifelse(is.null(phmmb2$fail_conv),phmmb2$sigma,NA),
                 ifelse(is.null(phmm$fail_conv),phmm$sigma,NA),
                 ifelse(is.null(phmm$fail_conv),phmm$sigma,NA)
               ),each=nrow(df)), 
      convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                    ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),
                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),
                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c( p$conv_problem,p2$conv_problem,  
                     pac$conv_problem,  pac2$conv_problem,   
                     ptva$conv_problem, ptva2$conv_problem, 
                     ptvb$conv_problem, ptvb2$conv_problem, 
                     ptvab$conv_problem,ptvab2$conv_problem,
                     phmma$conv_problem,phmma2$conv_problem,
                     phmmb$conv_problem,phmmb2$conv_problem,
                     phmm$conv_problem, phmm2$conv_problem   
                    ),each=nrow(df)))
    
    dfsig$pbias <- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias <- (dfsig$mode-dfsig$sim)

    
              
    #Smsy
    smsysim<-smsyCalc(dat$alpha,dat$beta)
  
    dfsmsy<- data.frame(parameter="smsy",
      iteration=u,
      scenario= simPars$scenario[a],
      method=rep(c(rep("MLE",16)),each=nrow(df)),
      model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
      prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
      by=rep(dat$year,16),
      sim=rep(smsysim,16),
      median=NA,
      mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$Smsy}, nrow(df)),
        rep(if(!is.null(p2$fail_conv)){NA}else{p2$Smsy}, nrow(df)),

        rep(if(!is.null(pac$fail_conv)){NA}else{pac$Smsy}, nrow(df)),
        rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$Smsy}, nrow(df)),

        if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smsy},
        if(!is.null(ptva2$fail_conv)){rep(NA, nrow(df))}else{ptva2$Smsy},

        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smsy},
        if(!is.null(ptvb2$fail_conv)){rep(NA, nrow(df))}else{ptvb2$Smsy},

        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smsy},
        if(!is.null(ptvab2$fail_conv)){rep(NA, nrow(df))}else{ptvab2$Smsy},

        if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$Smsy[phmma$regime]},
        if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$Smsy[phmma$regime]},

        if(!is.null(phmmb$fail_conv)){rep(NA, nrow(df))}else{phmmb$Smsy[phmmb$regime]},
        if(!is.null(phmmb2$fail_conv)){rep(NA, nrow(df))}else{phmmb2$Smsy[phmmb$regime]},

        if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$Smsy[phmm$regime]},
        if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$Smsy[phmm$regime]}

      ),    
        convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                        ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv)

                    ),each=nrow(df)),
        conv_warning=rep(c( p$conv_problem, p2$conv_problem,
                    pac$conv_problem, pac2$conv_problem,
                    ptva$conv_problem, ptva2$conv_problem,
                    ptvb$conv_problem, ptvb2$conv_problem,
                    ptvab$conv_problem, ptvab2$conv_problem,
                    phmma$conv_problem, phmma2$conv_problem,
                    phmmb$conv_problem, phmmb2$conv_problem,
                    phmm$conv_problem, phmm2$conv_problem
                    ),each=nrow(df))) 
  
  dfsmsy$pbias<- ((dfsmsy$mode-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (dfsmsy$mode-dfsmsy$sim)

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",16)),each=nrow(df)),
    model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
    prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
    by=rep(dat$year,16),
    sim=rep(unlist(mapply(sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),16),
    median=NA,
    mode=c(
      if(is.null(p$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="simple"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="simple"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="simple"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(p2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="simple"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="simple"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="simple"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},
 
      if(is.null(pac$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="autocorr"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="autocorr"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(pac2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="autocorr"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="autocorr"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},

      if(is.null(ptva$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwa"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(ptva2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwa"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},

      if(is.null(ptvb$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwb"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(ptvb2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwb"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},

      if(is.null(ptvab$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwab"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(ptvab2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="rwab"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},

      if(is.null(phmma$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(phmma2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmma"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmma"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="hmma"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},

      if(is.null(phmmb$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(phmmb2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmmb"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE],
           b=1/dfsmax$mode[dfsmax$model=="hmmb"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))},


      if(is.null(phmm$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"&dfa$prior==TRUE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"&dfsmsy$prior==TRUE], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"&dfsmax$prior==TRUE]))}else{rep(NA, nrow(df))},
      if(is.null(phmm2$fail_conv)){unlist(mapply(sGenCalc,
          a=dfa$mode[dfa$model=="hmmab"&dfa$method=="MLE"&dfa$prior==FALSE],
          Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"&dfsmsy$method=="MLE"&dfsmsy$prior==FALSE], 
          b=1/dfsmax$mode[dfsmax$model=="hmmab"&dfsmax$method=="MLE"&dfsmax$prior==FALSE]))}else{rep(NA, nrow(df))}
     ),
     
    convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                      ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv)
                    ),each=nrow(df)),
    conv_warning=rep(c( p$conv_problem, p2$conv_problem,
                     pac$conv_problem, pac2$conv_problem,
                     ptva$conv_problem, ptva2$conv_problem,
                     ptvb$conv_problem, ptvb2$conv_problem,
                     ptvab$conv_problem, ptvab2$conv_problem,
                     phmma$conv_problem, phmma2$conv_problem,
                     phmmb$conv_problem, phmmb2$conv_problem,
                     phmm$conv_problem, phmm2$conv_problem
                    ),each=nrow(df)))
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
 
  dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    scenario= simPars$scenario[a],
    method=rep(c(rep("MLE",16)),each=nrow(df)),
    model=rep(rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma","hmmb","hmmab"),each=2),each=nrow(df)),
    prior=rep(rep(c(TRUE,FALSE),each=nrow(df)),8),
    by=rep(dat$year,16),
    sim=rep(umsyCalc(dat$alpha),16),
    median=NA,
    mode=c(rep(if(!is.null(p$fail_conv)){NA}else{p$umsy}, nrow(df)),
           rep(if(!is.null(p2$fail_conv)){NA}else{p2$umsy}, nrow(df)),

           rep(if(!is.null(pac$fail_conv)){NA}else{pac$umsy}, nrow(df)),
           rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$umsy}, nrow(df)),

           if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$umsy},
           if(!is.null(ptva2$fail_conv)){rep(NA, nrow(df))}else{ptva2$umsy},

           if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$umsy},
           if(!is.null(ptvb2$fail_conv)){rep(NA, nrow(df))}else{ptvb2$umsy},

           if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$umsy},
           if(!is.null(ptvab2$fail_conv)){rep(NA, nrow(df))}else{ptvab2$umsy},

           if(!is.null(phmma$fail_conv)){rep(NA, nrow(df))}else{phmma$umsy[phmma$regime]},
           if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$umsy[phmma2$regime]},

           rep(if(!is.null(phmmb$fail_conv)){NA}else{phmmb$umsy}, nrow(df)),
           rep(if(!is.null(phmmb2$fail_conv)){NA}else{phmmb2$umsy}, nrow(df)),

           if(!is.null(phmm$fail_conv)){rep(NA, nrow(df))}else{phmm$umsy[phmm$regime]},
           if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$umsy[phmm2$regime]}

                  ), 
    convergence=rep(c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                      ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                    ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                    ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                    ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                    ifelse(is.null(phmm$fail_conv), phmm$model$convergence, phmm$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence, phmm2$fail_conv)

                    ),each=nrow(df)),
    conv_warning=rep(c( p$conv_problem, p2$conv_problem,
                    pac$conv_problem, pac2$conv_problem,
                    ptva$conv_problem, ptva2$conv_problem,
                    ptvb$conv_problem,  ptvb2$conv_problem,
                    ptvab$conv_problem, ptvab2$conv_problem,
                    phmma$conv_problem,  phmma2$conv_problem,
                    phmmb$conv_problem, phmmb2$conv_problem,
                    phmm$conv_problem, phmm2$conv_problem
                    ),each=nrow(df)))

     
    dfumsy$pbias<- ((dfumsy$mode-dfumsy$sim)/dfumsy$sim)*100
    dfumsy$bias<- (dfumsy$mode-dfumsy$sim)

    #AIC
    dfaic<- data.frame(parameter="AIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",16),
                       model=rep(c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab",
                                 "hmma","hmmb","hmmab"),each=2),
                       prior=rep(c(TRUE,FALSE),8),
                       by=NA,
                       sim=NA,
                       median=NA,
                       mode=c(ifelse(is.null(p$fail_conv),p$AICc, NA),
                              ifelse(is.null(p2$fail_conv),p2$AICc, NA),

                              ifelse(is.null(pac$fail_conv), pac$AICc, NA),
                              ifelse(is.null(pac2$fail_conv), pac2$AICc, NA),

                              ifelse(is.null(ptva$fail_conv), ptva$AICc, NA),
                              ifelse(is.null(ptva2$fail_conv), ptva2$AICc, NA),

                              ifelse(is.null(ptvb$fail_conv), ptvb$AICc, NA),
                              ifelse(is.null(ptvb2$fail_conv), ptvb2$AICc, NA),

                              ifelse(is.null(ptvab$fail_conv), ptvab$AICc, NA),
                              ifelse(is.null(ptvab2$fail_conv), ptvab2$AICc, NA),

                              ifelse(is.null(phmma$fail_conv), phmma$AICc, NA),
                              ifelse(is.null(phmma2$fail_conv), phmma2$AICc, NA),

                              ifelse(is.null(phmmb$fail_conv), phmmb$AICc, NA),
                              ifelse(is.null(phmmb2$fail_conv), phmmb2$AICc, NA),

                              ifelse(is.null(phmm$fail_conv), phmm$AICc, NA),
                              ifelse(is.null(phmm2$fail_conv), phmm2$AICc, NA)
                            ),
                       convergence=c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                              ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                              ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                              ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                              ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                              ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                              ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                              ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                              ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                              ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                              ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                              ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                              ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                              ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                              ifelse(is.null(phmm$fail_conv), phmm$model$convergence,phmm$fail_conv),
                              ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence, phmm2$fail_conv)
                               ),
                       conv_warning= c( p$conv_problem, p2$conv_problem,
                                      pac$conv_problem, pac2$conv_problem,
                                      ptva$conv_problem, ptva2$conv_problem,
                                      ptvb$conv_problem, ptvb2$conv_problem,
                                      ptvab$conv_problem, ptva2$conv_problem,
                                      phmma$conv_problem, phmma2$conv_problem,
                                      phmmb$conv_problem,  phmmb2$conv_problem,
                                      phmm$conv_problem, phmm2$conv_problem),
                       pbias=NA,
                       bias=NA)
    #BIC
    dfbic<- data.frame(parameter="BIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",16),
                       model=rep(c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab",
                                 "hmma","hmmb","hmmab"),each=2),
                       prior=rep(c(TRUE,FALSE),8),
                       by=NA,
                       sim=NA,
                       median=NA,
                       mode=c(ifelse(is.null(p$fail_conv),p$BIC, NA),
                              ifelse(is.null(p2$fail_conv),p2$BIC, NA),

                              ifelse(is.null(pac$fail_conv), pac$BIC, NA),
                              ifelse(is.null(pac2$fail_conv), pac2$BIC, NA),

                              ifelse(is.null(ptva$fail_conv), ptva$BIC, NA),
                              ifelse(is.null(ptva2$fail_conv), ptva2$BIC, NA),

                              ifelse(is.null(ptvb$fail_conv), ptvb$BIC, NA),
                              ifelse(is.null(ptvb2$fail_conv), ptvb2$BIC, NA),

                              ifelse(is.null(ptvab$fail_conv), ptvab$BIC, NA),
                              ifelse(is.null(ptvab2$fail_conv), ptvab2$BIC, NA),

                              ifelse(is.null(phmma$fail_conv), phmma$BIC, NA),
                               ifelse(is.null(phmma2$fail_conv), phmma2$BIC, NA),

                              ifelse(is.null(phmmb$fail_conv), phmmb$BIC, NA),
                              ifelse(is.null(phmmb2$fail_conv), phmmb2$BIC, NA),

                              ifelse(is.null(phmm$fail_conv), phmm$BIC, NA),
                              ifelse(is.null(phmm2$fail_conv), phmm2$BIC, NA)
                            ),
                       convergence=c(ifelse(is.null(p$fail_conv),p$model$convergence, p$fail_conv),
                                     ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),

                                     ifelse(is.null(pac$fail_conv), pac$model$convergence, pac$fail_conv),
                                     ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),

                                     ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                                     ifelse(is.null(ptva2$fail_conv), ptva2$model$convergence, ptva2$fail_conv),

                                     ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                                     ifelse(is.null(ptvb2$fail_conv), ptvb2$model$convergence, ptvb2$fail_conv),

                                     ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                                     ifelse(is.null(ptvab2$fail_conv), ptvab2$model$convergence, ptvab2$fail_conv),

                                     ifelse(is.null(phmma$fail_conv), phmma$model$convergence, phmma$fail_conv),
                                     ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),

                                     ifelse(is.null(phmmb$fail_conv), phmmb$model$convergence, phmmb$fail_conv),
                                     ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),

                                     ifelse(is.null(phmm$fail_conv), phmm$model$convergence, phmm$fail_conv),
                                     ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence, phmm2$fail_conv)
                                      ),
                       conv_warning= c( p$conv_problem,p2$conv_problem,
                                      pac$conv_problem, pac2$conv_problem,
                                      ptva$conv_problem, ptva2$conv_problem,
                                      ptvb$conv_problem, ptvb2$conv_problem,
                                      ptvab$conv_problem, ptvab2$conv_problem,
                                      phmma$conv_problem, phmma2$conv_problem,
                                      phmmb$conv_problem, phmmb2$conv_problem,
                                      phmm$conv_problem, phmm2$conv_problem),
                       pbias=NA,
                       bias=NA)

    
    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,
      dfaic,dfbic)

  return(dff)

}
