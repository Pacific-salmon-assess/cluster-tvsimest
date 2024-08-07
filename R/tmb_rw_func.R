
tmb_rw_func <- function(path=".",a, u) {
  
  dyn.load("src/srmodels/Ricker_tvlogb_centered.dll")
  dyn.load("src/srmodels/Ricker_tva_centered.dll")
  dyn.load("src/srmodels/Ricker_tva_tvb_centered.dll")

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
  
  dirpr<-matrix(c(2,1,1,2),2,2)

  
  ptva <- ricker_rw_TMB(data=df,tv.par='a',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
  ptvb <- ricker_rw_TMB(data=df, tv.par='b',sigb_p_sd=1,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
  ptvab <- ricker_rw_TMB(data=df, tv.par='both',sigb_p_sd=.4,logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
  

  #ptva$model$convergence
  #ptvb$model$convergence
  #ptvab$model$convergence

  #centered example
  initlm<-lm(logRS~S, data=df)

  tmb_paramsa <- list(alphao = max(initlm$coefficients[[1]],.5), 
                   logbeta=log(1/(max(df$S)*.5)),        
                   logsigobs = log(.6),
                   logsiga = log(.2),
                   adev=rep(0, length(df$S)-1))


  tmb_paramsb <- list(logbetao = log(1/(max(df$S)*.5)),
                   alpha   = max(initlm$coefficients[[1]],.5),                 
                   logsigobs = log(.6),
                   logsigb = log(.2),
                   bdev=rep(0, length(df$S)-1))


  tmb_paramsab <- list(logbetao = log(1/(max(df$S)*.5)),
                   alphao   = max(initlm$coefficients[[1]],.5),                 
                   logsigobs = log(.6),
                   logsiga = log(.5),
                   logsigb = log(.5),
                   bdev=rep(-1, length(df$S)-1),
                   adev=rep(0, length(df$S)-1))


  tmb_data_a <- list(
    obs_S = df$S,
    obs_logRS = df$logRS,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    siga_p_sd=1,
    logb_p_mean=logbeta_pr,
    logb_p_sd=logbeta_pr_sig
  )


  tmb_data_b <- list(
    obs_S = df$S,
    obs_logRS = df$logRS,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    sigb_p_sd=1,
    logb_p_mean=logbeta_pr,
    logb_p_sd=logbeta_pr_sig
  )


  tmb_data_ab <- list(
    obs_S = df$S,
    obs_logRS = df$logRS,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    sigb_p_sd=1,
    siga_p_sd=1,
    logb_p_mean=logbeta_pr,
    logb_p_sd=logbeta_pr_sig
  )


  obj_cent_a<-MakeADFun(tmb_data_a,tmb_paramsa,DLL="Ricker_tva_centered",random ="adev")
    newtonOption(obj_cent_a, smartsearch=FALSE)


  obj_cent_b<-MakeADFun(tmb_data_b,tmb_paramsb,DLL="Ricker_tvlogb_centered",random ="bdev")
    newtonOption(obj_cent_b, smartsearch=FALSE)


  obj_cent_ab<-MakeADFun(tmb_data_ab,tmb_paramsab,DLL="Ricker_tva_tvb_centered",random =c("bdev","adev"))
    newtonOption(obj_cent_ab, smartsearch=FALSE)

opt_cent_a<-nlminb(obj_cent_a$par,obj_cent_a$fn,obj_cent_a$gr)
opt_cent_b<-nlminb(obj_cent_b$par,obj_cent_b$fn,obj_cent_b$gr)
opt_cent_ab<-nlminb(obj_cent_ab$par,obj_cent_ab$fn,obj_cent_ab$gr)

rep_cent_a<-obj_cent_a$report()
rep_cent_b<-obj_cent_b$report()
rep_cent_ab<-obj_cent_ab$report()

sd_reporta <- TMB::sdreport(obj_cent_a)
conva <- get_convergence_diagnostics(sd_reporta)

sd_reportb <- TMB::sdreport(obj_cent_b)
convb <- get_convergence_diagnostics(sd_reportb)


sd_reportab <- TMB::sdreport(obj_cent_ab)
convab <- get_convergence_diagnostics(sd_reportab)


  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",8)),each=nrow(df)),
              model=rep(c(
                   "rwa","rwb","rwab",
                   "rwa_c","rwb_c","rwab_c"),each=nrow(df)),
              by=rep(dat$year,8),
              sim=rep(dat$alpha,8),
              median=NA,
              mode=c(
                    ptva$alpha,
                    rep(ptvb$alpha,nrow(df)),
                    ptvab$alpha,
                    rep_cent_a$alpha,
                    rep(rep_cent_b$alpha,nrow(df)),
                    rep_cent_ab$alpha,
                    ), 
              convergence=c(rep(c(
                    ptva$model$convergence,
                    ptvb$model$convergence,
                    ptvab$model$convergence, 
                    opt_cent_a$convergence,
                    opt_cent_b$convergence,
                    opt_cent_ab$convergence
                    ),each=nrow(df))),
              conv_warning=c(rep(c(
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                    conva$conv_problem,
                    convb$conv_problem,
                    convab$conv_problem
                    ),each=nrow(df))))
                    
  dfa$pbias <- ((dfa$mode-dfa$sim)/dfa$sim)*100
  dfa$bias <- (dfa$mode-dfa$sim)
  
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
      convergence=c(rep(c(p$model$convergence,
                    pac$model$convergence ,
                    ptva$model$convergence,
                    ptvb$model$convergence,
                    ptvab$model$convergence ,
                    phmma$model$convergence ,
                    phmmb$model$convergence ,
                    phmm$model$convergence
                    ),each=nrow(df))),
      conv_warning=c(rep(c( p$conv_problem,
                     pac$conv_problem,
                    ptva$conv_problem,
                    ptvb$conv_problem,
                    ptvab$conv_problem,
                     phmma$conv_problem,
                     phmmb$conv_problem,
                    phmm$conv_problem
                    ),each=nrow(df))))
      
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
        phmm$model$convergence + phmm$conv_problem),each=nrow(df)),
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
      mode=c(ptva$siga,ptvab$siga),
      convergence=c( ptva$model$convergence ,       
        ptvab$model$convergence ),
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
      mode=c(ptvb$sigb,ptvab$sigb),
      convergence=c( ptvb$model$convergence,        
        ptvab$model$convergence),
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
      mode=c(rep(p$Smsy,nrow(df)),
            rep(pac$Smsy,nrow(df)),
            ptva$Smsy,
            ptvb$Smsy,
            ptvab$Smsy,
            phmma$Smsy[phmma$regime],
            phmmb$Smsy[phmmb$regime],
            phmm$Smsy[phmm$regime]),    
        convergence=rep(c(p$model$convergence,
                    pac$model$convergence ,
                    ptva$model$convergence,
                    ptvb$model$convergence,
                    ptvab$model$convergence ,
                    phmma$model$convergence ,
                    phmmb$model$convergence ,
                    phmm$model$convergence
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
    convergence=rep(c(p$model$convergence,
                    pac$model$convergence ,
                    ptva$model$convergence,
                    ptvb$model$convergence,
                    ptvab$model$convergence ,
                    phmma$model$convergence ,
                    phmmb$model$convergence ,
                    phmm$model$convergence
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
    mode=c(rep(p$umsy, nrow(df)),
          rep(pac$umsy, nrow(df)),
           ptva$umsy,
          rep(ptvb$umsy, nrow(df)),
          ptvab$umsy,
          phmma$umsy[phmma$regime],
          rep(phmmb$umsy,nrow(df)),
          phmm$umsy[phmm$regime]),
    convergence=rep(c(p$model$convergence,
                    pac$model$convergence ,
                    ptva$model$convergence,
                    ptvb$model$convergence,
                    ptvab$model$convergence ,
                    phmma$model$convergence ,
                    phmmb$model$convergence ,
                    phmm$model$convergence
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
                       mode=c(p$AICc,
                             pac$AICc,
                             ptva$AICc,
                             ptvb$AICc,
                             ptvab$AICc,
                             phmma$AICc,
                             phmmb$AICc,
                             phmm$AICc),
                       convergence=c(p$model$convergence ,
                                     pac$model$convergence ,
                                     ptva$model$convergence ,
                                     ptvb$model$convergence,
                                     ptvab$model$convergence ,
                                     phmma$model$convergence ,
                                     phmmb$model$convergence ,
                                     phmm$model$convergence ),
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
                       mode=c(p$BIC,
                             pac$BIC,
                             ptva$BIC,
                             ptvb$BIC,
                             ptvab$BIC,
                             phmma$BIC,
                             phmmb$BIC,
                             phmm$BIC),
                       convergence=c(p$model$convergence ,
                                     pac$model$convergence ,
                                     ptva$model$convergence,
                                     ptvb$model$convergence,
                                     ptvab$model$convergence ,
                                     phmma$model$convergence ,
                                     phmmb$model$convergence,
                                     phmm$model$convergence),
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

    dff<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy,dfsiga,dfsigb,dfaic,dfbic,dflfo)

  return(dff)

}
