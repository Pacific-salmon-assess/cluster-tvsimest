
tmb_func_rw_comp <- function(path=".",a, u) {
  
  
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
 
  logbeta_pr_sig=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_pr=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
  

  ptva <-tryCatch({ ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,
                  logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001,useEDF=FALSE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvaEDF <-tryCatch({ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,
                  logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001,useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvac <-tryCatch({ricker_rw_TMB_centered(data=df,tv.par='a',sig_p_sd=1,
                    logb_p_mean=logbeta_pr, logb_p_sd=logbeta_pr_sig,   
                     deltaEDF=0.0001,useEDF=FALSE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvacEDF <-tryCatch({ricker_rw_TMB_centered(data=df,tv.par='a',sig_p_sd=1,
                    logb_p_mean=logbeta_pr, logb_p_sd=logbeta_pr_sig,   
                     deltaEDF=0.0001,useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvb <-tryCatch({ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001,useEDF=FALSE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvbEDF <-tryCatch({ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001,useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvbc <-tryCatch({ricker_rw_TMB_centered(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001,useEDF=FALSE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})
  
  ptvbcEDF <-tryCatch({ricker_rw_TMB_centered(data=df, tv.par='b',sigb_p_sd=1,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001,useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})
  
  ptvab <-tryCatch({ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvabEDF <-tryCatch({ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001, useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvabc <-tryCatch({ricker_rw_TMB_centered(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, deltaEDF=0.0001)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  ptvabcEDF <-tryCatch({ricker_rw_TMB_centered(data=df, tv.par='both',sigb_p_sd=.4,
                   logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig, 
                   deltaEDF=0.0001, useEDF=TRUE)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))})

  
  
  dfa<- data.frame(parameter="logalpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method="MLE",
              model=rep(rep(c("rwa","rwb","rwab"),each=2),each=nrow(df)),
              EDF=rep(rep(c(FALSE),6),each=nrow(df)),
              type=rep(rep(c("initial","centered"),3),each=nrow(df)),
              by=rep(dat$year,6),
              sim=rep(dat$alpha,6),
              median=NA,
              mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$alpha},
                    if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$logalpha},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{rep(ptvb$alpha, nrow(df))},
                    if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$logalpha},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$alpha},
                    if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$logalpha}), 
              convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
              conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,          
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
  #Smax
  dfsmax<- data.frame(parameter="smax",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(c(rep(c("rwa","rwb","rwab"),each=2)),each=nrow(df)),
      EDF=FALSE,
      type=rep(rep(c("initial","centered"),3),each=nrow(df)),
      by=rep(dat$year,6),
      sim=rep(1/dat$beta,6),
      median=NA,
      mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{rep(ptva$Smax, nrow(df))},
              if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$Smax},
              if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smax},
              if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$Smax},
              if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smax},
              if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$Smax}), 
      convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                     
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
      conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))
      
    dfsmax$pbias <- ((dfsmax$mode-dfsmax$sim)/dfsmax$sim)*100
    dfsmax$bias <- (dfsmax$mode-dfsmax$sim)
       
    #sigma
    dfsig<- data.frame(parameter="sigma",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(c(rep(c("rwa","rwb","rwab"),each=2)),each=nrow(df)),
      EDF=FALSE,
      type=rep(rep(c("initial","centered"),3),each=nrow(df)),
      by=rep(dat$year,6),
      sim=rep(dat$sigma,6),
      median=NA,
      mode=rep(c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$sig},
             if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$sig},
              if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$sig},
              if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$sig},
              if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$sig},
              if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$sig}),each=nrow(df)),
      convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
     conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))


    dfsig$pbias <- ((dfsig$mode-dfsig$sim)/dfsig$sim)*100
    dfsig$bias <- (dfsig$mode-dfsig$sim)

    #sigma_a
    
    dfsiga<- data.frame(parameter="sigma_a",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=c(rep(c("rwa","rwab"),each=2)),
      EDF=FALSE,
      type=rep(c("initial","centered"),2),
      by=NA,
      sim=NA,
      median=NA,
      mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$siga},
             if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$siga},
              if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$siga},
              if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$siga}),
      convergence=c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),
     conv_warning=c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ))


    dfsiga$pbias<- ((dfsiga$mode-dfsiga$sim)/dfsiga$sim)*100
    dfsiga$bias<- (dfsiga$mode-dfsiga$sim)

    #sigma_b
    dfsigb<- data.frame(parameter="sigma_b",
      iteration=u,
      scenario= simPars$scenario[a],
      method="MLE",
      model=rep(c("rwb","rwab"),each=2),
      EDF=FALSE,
      type=rep(c("initial","centered"),2),
      by=NA,
      sim=NA,
      median=NA,
      mode=c(if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$sigb},
              if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$sigb},
              if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$sigb},
              if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$sigb}),
      convergence=c(ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),
     conv_warning=c( ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ))

    dfsigb$pbias <- ((dfsigb$mode-dfsigb$sim)/dfsigb$sim)*100
    dfsigb$bias <- (dfsigb$mode-dfsigb$sim)
              
    #Smsy
    smsysim<-smsyCalc(dat$alpha,dat$beta)
  
    dfsmsy<- data.frame(parameter="smsy",
              iteration=u,
              scenario= simPars$scenario[a],
              method="MLE",
              model=rep(rep(c("rwa","rwb","rwab"),each=2),each=nrow(df)),
              EDF=rep(rep(c(FALSE),6),each=nrow(df)),
              type=rep(rep(c("initial","centered"),3),each=nrow(df)),
              by=rep(dat$year,6),
              sim=rep(smsysim,6),
              median=NA,
              mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$Smsy},
                    if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$Smsy},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{ptvb$Smsy},
                    if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$Smsy},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$Smsy},
                    if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$Smsy}), 
              convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
              conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,          
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))
                    
  dfsmsy$pbias<- ((dfsmsy$mode-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (dfsmsy$mode-dfsmsy$sim)

  
  #Sgen
  dfsgen <- data.frame(parameter="sgen",
              iteration=u,
              scenario= simPars$scenario[a],
              method="MLE",
              model=rep(rep(c("rwa","rwb","rwab"),each=2),each=nrow(df)),
              EDF=rep(rep(c(FALSE),6),each=nrow(df)),
              type=rep(rep(c("initial","centered"),3),each=nrow(df)),
              by=rep(dat$year,6),
              sim=rep(smsysim,6),
              median=NA,
              mode=c(if(is.null(ptva$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwa"&dfa$type=="initial"],
                      Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$type=="initial"], 
                      b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$type=="initial"]))}else{rep(NA, nrow(df))},
                    if(is.null(ptvac$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwa"&dfa$type=="centered"],
                      Smsy=dfsmsy$mode[dfsmsy$model=="rwa"&dfsmsy$type=="centered"], 
                      b=1/dfsmax$mode[dfsmax$model=="rwa"&dfsmax$type=="centered"]))}else{rep(NA, nrow(df))},
                    
                    if(is.null(ptvb$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwb"&dfa$type=="initial"],
                       Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$type=="initial"], 
                       b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$type=="initial"]))}else{rep(NA, nrow(df))},
                    if(is.null(ptvbc$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwb"&dfa$type=="centered"],
                       Smsy=dfsmsy$mode[dfsmsy$model=="rwb"&dfsmsy$type=="centered"], 
                       b=1/dfsmax$mode[dfsmax$model=="rwb"&dfsmax$type=="centered"]))}else{rep(NA, nrow(df))},

                    if(is.null(ptvab$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwab"&dfa$type=="initial"],
                      Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$type=="initial"], 
                      b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$type=="initial"]))}else{rep(NA, nrow(df))},
                    if(is.null(ptvabc$fail_conv)){unlist(mapply(sGenCalc,a=dfa$mode[dfa$model=="rwab"&dfa$type=="centered"],
                      Smsy=dfsmsy$mode[dfsmsy$model=="rwab"&dfsmsy$type=="centered"], 
                      b=1/dfsmax$mode[dfsmax$model=="rwab"&dfsmax$type=="centered"]))}else{rep(NA, nrow(df))}),
              convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
              conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,          
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))
                    
  
    dfsgen$pbias<- ((dfsgen$mode-dfsgen$sim)/dfsgen$sim)*100
    dfsgen$bias<- (dfsgen$mode-dfsgen$sim)
         
  #umsy
 
  dfumsy <- data.frame(parameter="umsy",
              iteration=u,
              scenario= simPars$scenario[a],
              method="MLE",
              model=rep(rep(c("rwa","rwb","rwab"),each=2),each=nrow(df)),
              EDF=rep(rep(c(FALSE),6),each=nrow(df)),
              type=rep(rep(c("initial","centered"),3),each=nrow(df)),
              by=rep(dat$year,6),
              sim=rep(dat$alpha,6),
              median=NA,
              mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$umsy},
                    if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$umsy},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{rep(ptvb$umsy, nrow(df))},
                    if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$umsy},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$umsy},
                    if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$umsy}), 
              convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
              conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,          
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)))

    dfumsy$pbias<- ((dfumsy$mode-dfumsy$sim)/dfumsy$sim)*100
    dfumsy$bias<- (dfumsy$mode-dfumsy$sim)


     #AIC
    dfaic<- data.frame(parameter="AIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",12),
                       model=rep(rep(c("rwa","rwb","rwab"),each=4),each=nrow(df)),
                       EDF=rep(rep(c(FALSE,TRUE),6),each=nrow(df)),
                       type=rep(rep(c("initial","initial","centered","centered"),3),each=nrow(df)),
                       by=rep(NA,12),
                       sim=rep(NA,12),
                       median=NA,
                       mode=c(ifelse(is.null(ptva$fail_conv), ptva$AICc, NA),
                              ifelse(is.null(ptva$fail_conv), ptvaEDF$AICc, NA),
                              ifelse(is.null(ptva$fail_conv), ptvac$AICc, NA),
                              ifelse(is.null(ptva$fail_conv), ptvacEDF$AICc, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvb$AICc, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbEDF$AICc, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbc$AICc, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbcEDF$AICc, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvab$AICc, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabEDF$AICc, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabc$AICc, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabcEDF$AICc, NA)),
                       convergence=c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                              ifelse(is.null(ptvaEDF$fail_conv), ptvaEDF$model$convergence, ptvaEDF$fail_conv),
                              ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                              ifelse(is.null(ptvacEDF$fail_conv), ptvacEDF$model$convergence, ptvacEDF$fail_conv),

                              ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                              ifelse(is.null(ptvbEDF$fail_conv), ptvbEDF$model$convergence, ptvbEDF$fail_conv),
                              ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                              ifelse(is.null(ptvbcEDF$fail_conv), ptvbcEDF$model$convergence, ptvbcEDF$fail_conv),

                              ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                              ifelse(is.null(ptvabEDF$fail_conv), ptvabEDF$model$convergence, ptvabEDF$fail_conv),
                              ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv),
                              ifelse(is.null(ptvabcEDF$fail_conv), ptvabcEDF$model$convergence, ptvabcEDF$fail_conv)),
              conv_warning=c( 
                    ptva$conv_problem,
                    ptvaEDF$conv_problem,
                    ptvac$conv_problem,  
                    ptvacEDF$conv_problem,           
                    ptvb$conv_problem,
                    ptvbEDF$conv_problem,
                    ptvbc$conv_problem,
                    ptvbcEDF$conv_problem,
                    ptvab$conv_problem,
                    ptvabEDF$conv_problem,
                    ptvabc$conv_problem,
                    ptvabcEDF$conv_problem),
                       pbias=rep(NA,12),
                       bias=rep(NA,12))
    #BIC
    dfbic<- data.frame(parameter="BIC",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",12),
                       model=rep(rep(c("rwa","rwb","rwab"),each=4),each=nrow(df)),
                       EDF=rep(rep(c(FALSE,TRUE),6),each=nrow(df)),
                       type=rep(rep(c("initial","initial","centered","centered"),3),each=nrow(df)),
                       by=rep(NA,12),
                       sim=rep(NA,12),
                       median=NA,
                        mode=c(ifelse(is.null(ptva$fail_conv), ptva$BIC, NA),
                              ifelse(is.null(ptva$fail_conv), ptvaEDF$BIC, NA),
                              ifelse(is.null(ptva$fail_conv), ptvac$BIC, NA),
                              ifelse(is.null(ptva$fail_conv), ptvacEDF$BIC, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvb$BIC, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbEDF$BIC, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbc$BIC, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbcEDF$BIC, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvab$BIC, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabEDF$BIC, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabc$BIC, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabcEDF$BIC, NA)),
                       convergence=c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                              ifelse(is.null(ptvaEDF$fail_conv), ptvaEDF$model$convergence, ptvaEDF$fail_conv),
                              ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                              ifelse(is.null(ptvacEDF$fail_conv), ptvacEDF$model$convergence, ptvacEDF$fail_conv),

                              ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                              ifelse(is.null(ptvbEDF$fail_conv), ptvbEDF$model$convergence, ptvbEDF$fail_conv),
                              ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                              ifelse(is.null(ptvbcEDF$fail_conv), ptvbcEDF$model$convergence, ptvbcEDF$fail_conv),

                              ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                              ifelse(is.null(ptvabEDF$fail_conv), ptvabEDF$model$convergence, ptvabEDF$fail_conv),
                              ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv),
                              ifelse(is.null(ptvabcEDF$fail_conv), ptvabcEDF$model$convergence, ptvabcEDF$fail_conv)),
                       conv_warning=c( 
                              ptva$conv_problem,
                              ptvaEDF$conv_problem,
                              ptvac$conv_problem,  
                              ptvacEDF$conv_problem,           
                              ptvb$conv_problem,
                              ptvbEDF$conv_problem,
                              ptvbc$conv_problem,
                              ptvbcEDF$conv_problem,
                              ptvab$conv_problem,
                              ptvabEDF$conv_problem,
                              ptvabc$conv_problem,
                              ptvabcEDF$conv_problem),
                       pbias=rep(NA,12),
                       bias=rep(NA,12))
                                 

     #EDF
    
    dfedf<- data.frame(parameter="EDF",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep("MLE",12),
                       model=rep(rep(c("rwa","rwb","rwab"),each=4),each=nrow(df)),
                       EDF=rep(rep(c(FALSE,TRUE),6),each=nrow(df)),
                       type=rep(rep(c("initial","initial","centered","centered"),3),each=nrow(df)),
                       by=rep(NA,12),
                       sim=rep(NA,12),
                       median=NA,
                        mode=c(ifelse(is.null(ptva$fail_conv), 4, NA),
                              ifelse(is.null(ptva$fail_conv), ptvaEDF$EDF, NA),
                              ifelse(is.null(ptva$fail_conv), 4, NA),
                              ifelse(is.null(ptva$fail_conv), ptvacEDF$EDF, NA),
                              ifelse(is.null(ptvb$fail_conv), 4, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbEDF$EDF, NA),
                              ifelse(is.null(ptvb$fail_conv), 4, NA),
                              ifelse(is.null(ptvb$fail_conv), ptvbcEDF$EDF, NA),
                              ifelse(is.null(ptvab$fail_conv), 5, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabEDF$EDF, NA),
                              ifelse(is.null(ptvab$fail_conv), 5, NA),
                              ifelse(is.null(ptvab$fail_conv), ptvabcEDF$EDF, NA)),
                       convergence=c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                              ifelse(is.null(ptvaEDF$fail_conv), ptvaEDF$model$convergence, ptvaEDF$fail_conv),
                              ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                              ifelse(is.null(ptvacEDF$fail_conv), ptvacEDF$model$convergence, ptvacEDF$fail_conv),

                              ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                              ifelse(is.null(ptvbEDF$fail_conv), ptvbEDF$model$convergence, ptvbEDF$fail_conv),
                              ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                              ifelse(is.null(ptvbcEDF$fail_conv), ptvbcEDF$model$convergence, ptvbcEDF$fail_conv),

                              ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                              ifelse(is.null(ptvabEDF$fail_conv), ptvabEDF$model$convergence, ptvabEDF$fail_conv),
                              ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv),
                              ifelse(is.null(ptvabcEDF$fail_conv), ptvabcEDF$model$convergence, ptvabcEDF$fail_conv)),
                       conv_warning=c( 
                              ptva$conv_problem,
                              ptvaEDF$conv_problem,
                              ptvac$conv_problem,  
                              ptvacEDF$conv_problem,           
                              ptvb$conv_problem,
                              ptvbEDF$conv_problem,
                              ptvbc$conv_problem,
                              ptvbcEDF$conv_problem,
                              ptvab$conv_problem,
                              ptvabEDF$conv_problem,
                              ptvabc$conv_problem,
                              ptvabcEDF$conv_problem),
                       pbias=rep(NA,12),
                       bias=rep(NA,12))


      dfgrad_i<- data.frame(parameter="grad_i",
              iteration=u,
              scenario= simPars$scenario[a],
              method="MLE",
              model=rep(rep(c("rwa","rwb","rwab"),each=2),each=nrow(df)),
              EDF=rep(rep(c(FALSE),6),each=nrow(df)),
              type=rep(rep(c("initial","centered"),3),each=nrow(df)),
              by=rep(dat$year,6),
              sim=rep(dat$alpha,6),
              median=NA,
              mode=c(if(!is.null(ptva$fail_conv)){rep(NA, nrow(df))}else{ptva$grad_i},
                    if(!is.null(ptvac$fail_conv)){rep(NA, nrow(df))}else{ptvac$umsy},
                    if(!is.null(ptvb$fail_conv)){rep(NA, nrow(df))}else{rep(ptvb$umsy, nrow(df))},
                    if(!is.null(ptvbc$fail_conv)){rep(NA, nrow(df))}else{ptvbc$umsy},
                    if(!is.null(ptvab$fail_conv)){rep(NA, nrow(df))}else{ptvab$umsy},
                    if(!is.null(ptvabc$fail_conv)){rep(NA, nrow(df))}else{ptvabc$umsy}), 
              convergence=rep(c(ifelse(is.null(ptva$fail_conv), ptva$model$convergence, ptva$fail_conv),
                    ifelse(is.null(ptvac$fail_conv), ptvac$model$convergence, ptvac$fail_conv),
                    ifelse(is.null(ptvb$fail_conv), ptvb$model$convergence, ptvb$fail_conv),
                    ifelse(is.null(ptvbc$fail_conv), ptvbc$model$convergence, ptvbc$fail_conv),
                    ifelse(is.null(ptvab$fail_conv), ptvab$model$convergence, ptvab$fail_conv),
                    ifelse(is.null(ptvabc$fail_conv), ptvabc$model$convergence, ptvabc$fail_conv)),each=nrow(df)),
              conv_warning=rep(c( 
                    ptva$conv_problem,
                    ptvac$conv_problem,          
                    ptvb$conv_problem,
                    ptvbc$conv_problem,
                    ptvab$conv_problem,
                    ptvabc$conv_problem
                    ),each=nrow(df)),
              pbias=rep(NA,12),
              bias=rep(NA,12))
    
  

  
  

    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfsiga,dfsigb,
      dfaic,dfbic,dfedf,dfgrad_i)

    print(paste("a is ",a," u is ", u))

  return(dff)

}
