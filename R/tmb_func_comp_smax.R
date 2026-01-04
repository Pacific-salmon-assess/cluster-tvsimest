
tmb_func_comp_smax <- function(path=".",a, u) {
  
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


  p3 <- tryCatch({ ricker_TMB(data=df,Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  p2 <- tryCatch({ ricker_TMB_deprecated(data=df,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  pac3 <- tryCatch({ricker_TMB(data=df, AC=TRUE,Smax_mean=Smax_mean,
                 Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  pac2 <- tryCatch({ricker_TMB_deprecated(data=df, AC=TRUE,logb_p_mean=logbeta_pr,
                 logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptva <- tryCatch({ricker_rw_TMB_logb(data=df,tv.par='a',logb_p_mean=logbeta_pr,
                  logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptva3 <- tryCatch({ricker_rw_TMB(data=df,tv.par='a',Smax_mean=Smax_mean,
                 Smax_sd=Smax_sd, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvb <- tryCatch({ricker_rw_TMB_logb(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvb3 <- tryCatch({ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
                   Smax_mean=Smax_mean, Smax_sd=Smax_sd, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})
  
  ptvab <- tryCatch({ricker_rw_TMB_logb(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvab3 <- tryCatch({ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,
                   Smax_mean=Smax_mean, Smax_sd=Smax_sd, deltaEDF=0.0001, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmma3 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='a', dirichlet_prior=dirpr,
                     Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})


  phmmb3 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='b', dirichlet_prior=dirpr,
                     Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmm3 <- tryCatch({ricker_hmm_TMB2(data=df, tv.par='both', dirichlet_prior=dirpr,
                   Smax_mean=Smax_mean,Smax_sd=Smax_sd, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )
  

  phmma2 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='a', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmmb2 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='b', dirichlet_prior=dirpr,
                    logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  phmm2 <- tryCatch({ricker_hmm_TMB2_logb(data=df, tv.par='both', dirichlet_prior=dirpr,
                  logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, silent=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))} )

  dfa <- data.frame(parameter="logalpha",
              iteration=u,
              scenario= simPars$scenario[a],
              version=rep(c(rep(c("logb","smax"),16)),each=nrow(df)),
              model=rep(c("simple","simple",
                   "autocorr","autocorr",
                   "rwa","rwa",
                   "rwb","rwb",
                   "rwab","rwab",
                   "hmma","hmma",
                   "hmmb","hmmb",
                   "hmmab","hmmab"),each=nrow(df)),
              by=rep(dat$year,16),
              sim=rep(dat$alpha,16),
              median=NA,
              mode=c(rep(if(!is.null(p2$fail_conv)){NA}else{p2$logalpha}, nrow(df)),
                    rep(if(!is.null(p3$fail_conv)){NA}else{p3$logalpha}, nrow(df)),
                    rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$logalpha}, nrow(df)),
                    rep(if(!is.null(pac3$fail_conv)){NA}else{pac3$logalpha}, nrow(df)),
                    if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$logalpha},
                    if(!is.null(ptva3$fail_conv)){rep(NA, nrow(df))}else{ptva3$logalpha},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$logalpha},
                    if(!is.null(ptvb3$fail_conv)){rep(NA, nrow(df))}else{ptvb3$logalpha},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$logalpha},
                    if(!is.null(ptvab3$fail_conv)){rep(NA, nrow(df))}else{ptvab3$logalpha},
                   if(!is.null(phmma2$fail_conv)){rep(NA, nrow(df))}else{phmma2$logalpha[phmma2$regime]},
                   if(!is.null(phmma3$fail_conv)){rep(NA, nrow(df))}else{phmma3$logalpha[phmma3$regime]},
                    rep(if(!is.null(phmmb2$fail_conv)){NA}else{phmmb2$logalpha}, nrow(df)),
                    rep(if(!is.null(phmmb3$fail_conv)){NA}else{phmmb3$logalpha}, nrow(df)),
                    if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$logalpha[phmm2$regime]},
                    if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$logalpha[phmm3$regime]}
                  ), 
              convergence=rep(c(
                    ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),
                    ifelse(is.null(p3$fail_conv),p3$model$convergence, p3$fail_conv),
                    
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),
                    ifelse(is.null(pac3$fail_conv), pac3$model$convergence, pac3$fail_conv),
                    
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva3$fail_conv), ptva3$model$convergence, ptva3$fail_conv),
                    
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb3$fail_conv), ptvb3$model$convergence, ptvb3$fail_conv),
                    
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab3$fail_conv), ptvab3$model$convergence, ptvab3$fail_conv),
 
                    
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
              conv_warning=rep(c(p2$conv_problem,
                p3$conv_problem,
                    pac2$conv_problem,
                    pac3$conv_problem,
                    ptva$conv_problem,
                    ptva3$conv_problem,
                    ptvb$conv_problem,
                    ptvb3$conv_problem,
                    ptvab$conv_problem,
                    ptvab3$conv_problem,
                    phmma2$conv_problem,
                    phmma3$conv_problem,
                    phmmb2$conv_problem,
                    phmmb3$conv_problem,
                    phmm2$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax <- data.frame(parameter="Smax",
      iteration=u,
      scenario= simPars$scenario[a],
      version=rep(c(rep(c("logb","smax"),16)),each=nrow(df)),
              model=rep(c("simple","simple",
                   "autocorr","autocorr",
                   "rwa","rwa",
                   "rwb","rwb",
                   "rwab","rwab",
                   "hmma","hmma",
                   "hmmb","hmmb",
                   "hmmab","hmmab"),each=nrow(df)),
      by=rep(dat$year,16),
      sim=rep(1/dat$beta,16),
      median=NA,
      mode=c(
       rep(if(!is.null(p2$fail_conv)){NA}else{p2$Smax}, nrow(df)),
       rep(if(!is.null(p3$fail_conv)){NA}else{p3$Smax}, nrow(df)),
        rep(if(!is.null(pac2$fail_conv)){NA}else{pac2$Smax}, nrow(df)),
         rep(if(!is.null(pac3$fail_conv)){NA}else{pac3$Smax}, nrow(df)),
        if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smax},
        if(!is.null(ptva3$fail_conv)){rep(NA, nrow(df))}else{rep(ptva3$Smax,nrow(df))},
        if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smax},
         if(!is.null(ptvb3$fail_conv)){rep(NA, nrow(df))}else{ptvb3$Smax},
        if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smax},
        if(!is.null(ptvab3$fail_conv)){rep(NA, nrow(df))}else{ptvab3$Smax},
        rep(if(!is.null(phmma2$fail_conv)){NA}else{phmma2$Smax}, nrow(df)),
        rep(if(!is.null(phmma3$fail_conv)){NA}else{phmma3$Smax}, nrow(df)),
        if(!is.null(phmmb3$fail_conv)){rep(NA, nrow(df))}else{phmmb3$Smax[phmmb3$regime]},
           if(!is.null(phmmb2$fail_conv)){rep(NA, nrow(df))}else{phmmb2$Smax[phmmb2$regime]},
     if(!is.null(phmm2$fail_conv)){rep(NA, nrow(df))}else{phmm2$Smax[phmm2$regime]},
        if(!is.null(phmm3$fail_conv)){rep(NA, nrow(df))}else{phmm3$Smax[phmm3$regime]}
      ),
      convergence=rep(c(
                    ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),
                    ifelse(is.null(p3$fail_conv),p3$model$convergence, p3$fail_conv),
                    
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),
                    ifelse(is.null(pac3$fail_conv), pac3$model$convergence, pac3$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva3$fail_conv), ptva3$model$convergence, ptva3$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb3$fail_conv), ptvb3$model$convergence, ptvb3$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab3$fail_conv), ptvab3$model$convergence, ptvab3$fail_conv),
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c(
                     p2$conv_problem,
                      p3$conv_problem,
                     pac2$conv_problem,
                      pac3$conv_problem,
                    ptva$conv_problem,
                    ptva3$conv_problem,
                    ptvb$conv_problem,
                    ptvb3$conv_problem,
                    ptvab$conv_problem,
                    ptvab3$conv_problem,
                    phmma2$conv_problem,
                    phmma3$conv_problem,
                    phmmb2$conv_problem,
                    phmmb3$conv_problem,
                    phmm2$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      version=rep(c(rep(c("logb","smax"),16)),each=nrow(df)),
      model=rep(c( "simple","simple",
                   "autocorr","autocorr",
                   "rwa","rwa",
                   "rwb","rwb",
                   "rwab","rwab",
                   "hmma","hmma",
                   "hmmb","hmmb",
                   "hmmab","hmmab"),each=nrow(df)),
      by=rep(dat$year,16),
      sim=rep(dat$sigma,16),
      median=NA,
      mode=rep(c(ifelse(is.null(p2$fail_conv),p2$sig,NA),
        ifelse(is.null(p3$fail_conv),p3$sigma,NA),
                 ifelse(is.null(pac2$fail_conv),pac2$sig,NA),
                 ifelse(is.null(pac3$fail_conv),pac3$sigma,NA),
                 ifelse(is.null(ptva$fail_conv),ptva$sigma,NA),
                 ifelse(is.null(ptva3$fail_conv),ptva$sigma,NA),
                 ifelse(is.null(ptvb$fail_conv),ptvb$sigma,NA),
                 ifelse(is.null(ptvb3$fail_conv),ptvb$sigma,NA),
                 ifelse(is.null(ptvab$fail_conv),ptvab$sigma,NA),

                 ifelse(is.null(ptvab3$fail_conv),ptvab$sigma,NA),

        
                 ifelse(is.null(phmma2$fail_conv),phmma2$sigma,NA),
                 ifelse(is.null(phmma3$fail_conv),phmma3$sigma,NA),
                 ifelse(is.null(phmmb2$fail_conv),phmmb2$sigma,NA),
                 ifelse(is.null(phmmb3$fail_conv),phmmb3$sigma,NA),
                 ifelse(is.null(phmm2$fail_conv),phmm2$sigma,NA),
                 ifelse(is.null(phmm3$fail_conv),phmm3$sigma,NA)
               ),each=nrow(df)), 
      convergence=rep(c(ifelse(is.null(p2$fail_conv),p2$model$convergence, p2$fail_conv),
        ifelse(is.null(p3$fail_conv),p3$model$convergence, p3$fail_conv),
                    ifelse(is.null(pac2$fail_conv), pac2$model$convergence, pac2$fail_conv),
                    ifelse(is.null(pac3$fail_conv), pac3$model$convergence, pac3$fail_conv),
                    ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptva3$fail_conv), ptva3$model$convergence, ptva3$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvb3$fail_conv), ptvb3$model$convergence, ptvb3$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvab3$fail_conv), ptvab3$model$convergence, ptvab3$fail_conv),
        
                    ifelse(is.null(phmma2$fail_conv), phmma2$model$convergence, phmma2$fail_conv),
                    ifelse(is.null(phmma3$fail_conv), phmma3$model$convergence, phmma3$fail_conv),
                    ifelse(is.null(phmmb2$fail_conv), phmmb2$model$convergence, phmmb2$fail_conv),
                    ifelse(is.null(phmmb3$fail_conv), phmmb3$model$convergence, phmmb3$fail_conv),
                    ifelse(is.null(phmm2$fail_conv), phmm2$model$convergence,phmm2$fail_conv),
                    ifelse(is.null(phmm3$fail_conv), phmm3$model$convergence,phmm3$fail_conv)
                    ),each=nrow(df)),
      conv_warning=rep(c(p3$conv_problem,
        p3$conv_problem,
                     pac2$conv_problem,
                      pac3$conv_problem,
                    ptva$conv_problem,
                    ptva3$conv_problem,
                    ptvb$conv_problem,
                    ptvb3$conv_problem,
                    ptvab$conv_problem,
                    ptvab3$conv_problem,
                     phmma2$conv_problem,
                     phmmb2$conv_problem,
                    phmm2$conv_problem,
                    phmma3$conv_problem,
                     phmmb3$conv_problem,
                    phmm3$conv_problem
                    ),each=nrow(df)))
    
    dfsig$pbias <- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias <- (dfsig$mode-dfsig$sim)

              
   

    dff<-rbind(dfa,dfsmax,dfsig)

  return(dff)

}
