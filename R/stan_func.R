stan_func<- function(path=".", a,u){
  
  
  simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
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
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=max(dat$obsSpawners)*.5,
             pSmax_sig=max(dat$obsSpawners)*.5,
             psig_b=max(dat$obsSpawners)*.5
  )

 

  
  #create folder to hold temp files
#  dir.create(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u,sep=''))
  
#  ls=list.files("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/")
#  if(length(ls)>5){
 #   unlink(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u-2,
#                 '/*',sep=''))
#}
  
  options(mc.cores = 5)
  print(paste("a is", a))
  print(paste("u is", u))
  
  #
  print("simple")
  f1 <- mod1$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f1_ip<-f1$summary()
  conv_f1_ip <- check_stan_conv(stansum=f1_ip)



  
   print("autocorr")
  f2 <- mod2$sample(data=df,
                    seed = 123,
                   chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f2_ip<-f2$summary()
  conv_f2_ip <- check_stan_conv(stansum=f2_ip)
  
  print("rwa")
  f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f3_ip<-f3$summary()
  conv_f3_ip <- check_stan_conv(stansum=f3_ip)
  
  print("rwb")
  f4 <- mod4$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f4_ip<-f4$summary()
  conv_f4_ip <- check_stan_conv(stansum=f4_ip)
  
  print("rwab")
  f5 <- mod5$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f5_ip<-f5$summary()
  conv_f5_ip <- check_stan_conv(stansum=f5_ip)
  print("hmma")
  f6 <- mod6$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f6_ip <- f6$summary()
  conv_f6_ip <- check_stan_conv(stansum=f6_ip)
  
  print("hmmb")
  f7 <- mod7$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f7_ip <- f7$summary()
  conv_f7_ip <- check_stan_conv(stansum=f7_ip)

  print("hmmab")
  f8 <- mod8$sample(data=df,
                    seed = 123,
                   chains = 6, 
                    iter_warmup = 2000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  f8_ip <- f8$summary()
  conv_f8_ip <- check_stan_conv(stansum=f8_ip)
  
  
  
  #Max. prod
  
  
  drf1<-as.data.frame(f1$draws( format="matrix"))
  drf2<-as.data.frame(f2$draws( format="matrix"))
  drf3<-as.data.frame(f3$draws( format="matrix"))
  drf4<-as.data.frame(f4$draws( format="matrix"))
  drf5<-as.data.frame(f5$draws( format="matrix"))
  drf6<-as.data.frame(f6$draws( format="matrix"))
  drf7<-as.data.frame(f7$draws( format="matrix"))
  drf8<-as.data.frame(f8$draws( format="matrix"))

  
  loga_f1 <- postmode(x=drf1$log_a)
  loga_f2 <- postmode(x=drf2$log_a)
  loga_f3 <- apply(drf3[,grep("log_a\\[",colnames(drf3))],2,postmode)
  loga_f4 <- postmode(x=drf4$log_a)
  loga_f5 <- apply(drf5[,grep("log_a\\[",colnames(drf5))],2,postmode)
  
  loga_f6_regime <- apply(drf6[,grep("log_a\\[",colnames(drf6))],2,postmode)
  zstar_f6 <- apply(drf6[,grep("zstar\\[",colnames(drf6))],2,function(x)which.max(tabulate(as.integer(x))))
  loga_f6 <- loga_f6_regime[zstar_f6]
  loga_f6_conv <- 
    conv_f6_ip[[2]]$sumconv[grep('log_a\\[',conv_f6_ip[[2]]$variable)][f6_ip$median[grep('zstar\\[',f6_ip$variable)]] +
    conv_f6_ip[[2]]$sumconv[grep('zstar\\[',conv_f6_ip[[2]]$variable)]


  loga_f7 <- postmode(x=drf7$log_a)

  loga_f8_regime <- apply(drf8[,grep("log_a\\[",colnames(drf8))],2,postmode)
  zstar_f8 <- apply(drf8[,grep("zstar\\[",colnames(drf8))],2,function(x)which.max(tabulate(as.integer(x))))
  loga_f8 <- loga_f8_regime[zstar_f8] 
  loga_f8_conv <- 
    conv_f8_ip[[2]]$sumconv[grep('log_a\\[',conv_f8_ip[[2]]$variable)][f8_ip$median[grep('zstar\\[',f8_ip$variable)]] +
    conv_f8_ip[[2]]$sumconv[grep('zstar\\[',conv_f8_ip[[2]]$variable)]

  phmma_alpha_regime_median<-f6_ip$median[grep('log_a\\[',f6_ip$variable)]
  phmma_zstar_regime_median<-f6_ip$median[grep('zstar\\[',f6_ip$variable)]
  phmma_alpha_median <-phmma_alpha_regime_median[phmma_zstar_regime_median]

  phmmab_alpha_regime_median<- f8_ip$median[grep('log_a\\[',f8_ip$variable)]
  phmmab_zstar_regime_median<- f8_ip$median[grep('zstar\\[',f8_ip$variable)]
  phmmab_alpha_median <- phmmab_alpha_regime_median[phmmab_zstar_regime_median]
  #plot(dat$year,phmmab_alpha_median)

  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                   model=rep(c("simple",
                               "autocorr",
                               "rwa","rwb","rwab",
                               "hmma","hmmb","hmmab"),each=nrow(dat)),
                   by=rep(dat$year,8),
                   sim=rep(dat$alpha,8),
                   median=c(rep(f1_ip$median[f1_ip$variable=='log_a'],nrow(dat)),
                         rep(f2_ip$median[f2_ip$variable=='log_a'],nrow(dat)),
                         f3_ip$median[grep('log_a\\[',f3_ip$variable)],
                         rep(f4_ip$median[f4_ip$variable=='log_a'],nrow(dat)),
                         f5_ip$median[grep('log_a\\[',f5_ip$variable)],
                         phmma_alpha_median,
                         rep(f7_ip$median[f7_ip$variable=='log_a'],nrow(dat)),
                         phmmab_alpha_median),
                   mode=c(rep(loga_f1,nrow(dat)),
                         rep(loga_f2,nrow(dat)),
                         loga_f3,
                         rep(loga_f4,nrow(dat)),
                         loga_f5,
                         loga_f6,
                         rep(loga_f7,nrow(dat)),
                         loga_f8),
                   #need to get in the harsher convergence criteria
                   convergence= c(rep(conv_f1_ip[[2]]$sumconv[conv_f1_ip[[2]]$variable=="log_a"],nrow(dat)),
                        rep(conv_f2_ip[[2]]$sumconv[conv_f2_ip[[2]]$variable=="log_a"],nrow(dat)),
                        conv_f3_ip[[2]]$sumconv[grep('log_a\\[',conv_f3_ip[[2]]$variable)],
                        rep(conv_f4_ip[[2]]$sumconv[conv_f4_ip[[2]]$variable=="log_a"],nrow(dat)),
                        conv_f5_ip[[2]]$sumconv[grep('log_a\\[',conv_f5_ip[[2]]$variable)],
                        loga_f6_conv,
                        rep(conv_f7_ip[[2]]$sumconv[conv_f7_ip[[2]]$variable=="log_a"],nrow(dat)),
                        loga_f8_conv),
                   conv_warning=NA)


  dfa$pbias<- ((as.numeric(dfa$mode)-dfa$sim)/dfa$sim)*100
  dfa$bias<- (as.numeric(dfa$mode)-dfa$sim)
 

  #Smax
  smax_f1 <- postmode(x=drf1$S_max)
  smax_f2 <- postmode(x=drf2$S_max)
  smax_f3 <- postmode(x=drf3$S_max)

  smax_f4 <- apply(drf4[,grep("Smax\\[",colnames(drf4))],2,postmode)
  smax_f5 <- apply(drf5[,grep("Smax\\[",colnames(drf5))],2,postmode)

  smax_f6 <- postmode(x=drf6$S_max)

  smax_f7_regime <- apply(drf7[,grep("S_max\\[",colnames(drf7))],2,postmode)
  zstar_f7 <- apply(drf7[,grep("zstar\\[",colnames(drf7))],2,function(x)which.max(tabulate(as.integer(x))))
  smax_f7 <- smax_f7_regime[zstar_f7]

  smax_f7_conv <- 
    conv_f7_ip[[2]]$sumconv[grep('S_max\\[',conv_f7_ip[[2]]$variable)][f7_ip$median[grep('zstar\\[',f7_ip$variable)]] +
    conv_f7_ip[[2]]$sumconv[grep('zstar\\[',conv_f7_ip[[2]]$variable)]


  smax_f8_regime <- apply(drf8[,grep("S_max\\[",colnames(drf8))],2,postmode)
  smax_f8 <- smax_f8_regime[zstar_f8]
 
  smax_f8_conv <- 
    conv_f8_ip[[2]]$sumconv[grep('S_max\\[',conv_f8_ip[[2]]$variable)][f8_ip$median[grep('zstar\\[',f8_ip$variable)]] +
    conv_f8_ip[[2]]$sumconv[grep('zstar\\[',conv_f8_ip[[2]]$variable)]

  
  phmmb_smax_regime_median<-f7_ip$median[grep('S_max\\[',f7_ip$variable)]
  phmmb_zstar_regime_median<-f7_ip$median[grep('zstar\\[',f7_ip$variable)]
  phmmb_smax_median <-phmmb_smax_regime_median[phmmb_zstar_regime_median]

  phmmab_smax_regime_median<- f8_ip$median[grep('S_max\\[',f8_ip$variable)]
  phmmab_smax_median <- phmmab_smax_regime_median[phmmab_zstar_regime_median]


  dfsmax<- data.frame(parameter="smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma","hmmb","hmmab"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(1/dat$beta,8),
                      median=c(rep(f1_ip$median[f1_ip$variable=='S_max'],nrow(dat)),
                        rep(f2_ip$median[f2_ip$variable=='S_max'],nrow(dat)),
                        rep(f3_ip$median[f3_ip$variable=='S_max'],nrow(dat)),
                        f4_ip$median[grep('Smax\\[',f4_ip$variable)],
                        f5_ip$median[grep('Smax\\[',f5_ip$variable)],
                        rep(f6_ip$median[f6_ip$variable=='S_max'],nrow(dat)),
                        phmmb_smax_median,
                        phmmab_smax_median),
                      mode=c(rep(smax_f1,nrow(dat)),
                         rep(smax_f2,nrow(dat)),
                         rep(smax_f3,nrow(dat)),
                         smax_f4,
                         smax_f5,
                         rep(smax_f6,nrow(dat)),
                         smax_f7,
                         smax_f8),
                      convergence=
                      c(rep(conv_f1_ip[[2]]$sumconv[conv_f1_ip[[2]]$variable=="S_max"],nrow(dat)),
                        rep(conv_f2_ip[[2]]$sumconv[conv_f2_ip[[2]]$variable=="S_max"],nrow(dat)),
                        rep(conv_f3_ip[[2]]$sumconv[conv_f3_ip[[2]]$variable=="S_max"],nrow(dat)),
                        conv_f4_ip[[2]]$sumconv[grep('Smax\\[',conv_f4_ip[[2]]$variable)],
                        conv_f5_ip[[2]]$sumconv[grep('Smax\\[',conv_f5_ip[[2]]$variable)],
                        rep(conv_f6_ip[[2]]$sumconv[conv_f6_ip[[2]]$variable=="S_max"],nrow(dat)),
                        smax_f7_conv,
                        smax_f8_conv),
                      conv_warning=NA)
                      
  dfsmax$pbias <- ((as.numeric(dfsmax$mode)-dfsmax$sim)/dfsmax$sim)*100
  dfsmax$bias <- (as.numeric(dfsmax$mode)-dfsmax$sim)


  #sigma -obs error
  sig_f1 <- postmode(x=drf1$sigma)
  sig_f2 <- postmode(x=drf2$sigma)
  sig_f3 <- postmode(x=drf3$sigma)
  sig_f4 <- postmode(x=drf4$sigma)
  sig_f5 <- postmode(x=drf5$sigma)
  sig_f6 <- postmode(x=drf6$sigma)
  sig_f7 <- postmode(x=drf7$sigma)
  sig_f8 <- postmode(x=drf8$sigma)

 
 


  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method="MCMC",
                     model=rep(c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma","hmmb","hmmab"),each=nrow(dat)),
                     by=rep(dat$year,8),
                     sim=rep(dat$sigma,8),
                     median=rep(c(f1_ip$median[f1_ip$variable=='sigma'],
                      f2_ip$median[f2_ip$variable=='sigma'],
                      f3_ip$median[f3_ip$variable=='sigma'],
                      f4_ip$median[f4_ip$variable=='sigma'],
                      f5_ip$median[f5_ip$variable=='sigma'],
                      f6_ip$median[f6_ip$variable=='sigma'],
                      f7_ip$median[f7_ip$variable=='sigma'],
                      f8_ip$median[f8_ip$variable=='sigma']),each=nrow(dat)),
                     mode=rep(c(sig_f1,
                      sig_f2,
                      sig_f3,
                      sig_f4,
                      sig_f5,
                      sig_f6,
                      sig_f7,
                      sig_f8),each=nrow(dat)),
                    convergence= rep(c(conv_f1_ip[[2]]$sumconv[conv_f1_ip[[2]]$variable=='sigma'],
                      conv_f2_ip[[2]]$sumconv[conv_f2_ip[[2]]$variable=='sigma'],
                      conv_f3_ip[[2]]$sumconv[conv_f3_ip[[2]]$variable=='sigma'],
                      conv_f4_ip[[2]]$sumconv[conv_f4_ip[[2]]$variable=='sigma'],
                      conv_f5_ip[[2]]$sumconv[conv_f5_ip[[2]]$variable=='sigma'],
                      conv_f6_ip[[2]]$sumconv[conv_f6_ip[[2]]$variable=='sigma'],
                      conv_f7_ip[[2]]$sumconv[conv_f7_ip[[2]]$variable=='sigma'],
                      conv_f8_ip[[2]]$sumconv[conv_f8_ip[[2]]$variable=='sigma']),each=nrow(dat)),
                    conv_warning=NA)
  
  dfsig$pbias<- ((as.numeric(dfsig$mode)-dfsig$sim)/dfsig$sim)*100
  dfsig$bias<- (as.numeric(dfsig$mode)-dfsig$sim)
  
  #sigma a

  siga_f3 <- postmode(x=drf3$"sigma_a")
  siga_f5 <- postmode(x=drf5$"sigma_a")
 
  dfsiga<- data.frame(parameter="sigma_a",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwa","rwab"),
                      by=NA,
                      sim=NA,
                      median=c(f3_ip$median[f3_ip$variable=='sigma_a'],
                            f5_ip$median[f5_ip$variable=='sigma_a']),
                      mode=c(siga_f3,
                        siga_f5),
                      convergence=c(conv_f3_ip[[2]]$sumconv[conv_f3_ip[[2]]$variable=='sigma_a'],
                        conv_f5_ip[[2]]$sumconv[conv_f5_ip[[2]]$variable=='sigma_a']),
                      conv_warning=NA
                        )

  dfsiga$pbias<- NA
  dfsiga$bias<- NA
  
  #sigma b
  sigb_f4 <- postmode(x=drf4$"sigma_b")
  sigb_f5 <- postmode(x=drf5$"sigma_b")
 
  dfsigb<- data.frame(parameter="sigma_b",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwb","rwab"),
                      by=NA,
                      sim=NA,
                      median=c(f4_ip$median[f4_ip$variable=='sigma_b'],
                        f5_ip$median[f5_ip$variable=='sigma_b']),
                      mode=c(sigb_f4,
                          sigb_f5),
                      convergence=c(conv_f4_ip[[2]]$sumconv[conv_f4_ip[[2]]$variable=='sigma_b'],
                        conv_f5_ip[[2]]$sumconv[conv_f5_ip[[2]]$variable=='sigma_b']),
                      conv_warning=NA)

                     
  dfsigb$pbias<- NA
  dfsigb$bias<- NA
  

  #S_msy
  Smsy_f1 <- postmode(x=drf1$"S_msy")
  Smsy_f2 <- postmode(x=drf2$"S_msy")

  Smsy_f3 <- apply(drf3[,grep("S_msy\\[",colnames(drf3))],2,postmode)
  Smsy_f4 <- apply(drf4[,grep("S_msy\\[",colnames(drf4))],2,postmode)
  Smsy_f5 <- apply(drf5[,grep("S_msy\\[",colnames(drf5))],2,postmode)
  

  Smsy_f6_regime <- apply(drf6[,grep("S_msy\\[",colnames(drf6))],2,postmode)
  Smsy_f6 <- Smsy_f6_regime[zstar_f6]
  Smsy_f6_conv <- 
    conv_f6_ip[[2]]$sumconv[grep('S_msy\\[',conv_f6_ip[[2]]$variable)][f6_ip$median[grep('zstar\\[',f6_ip$variable)]] +
    conv_f6_ip[[2]]$sumconv[grep('zstar\\[',conv_f6_ip[[2]]$variable)]

  Smsy_f7_regime <- apply(drf7[,grep("S_msy\\[",colnames(drf7))],2,postmode)
  Smsy_f7 <- Smsy_f7_regime[zstar_f7]
  Smsy_f7_conv <- 
    conv_f7_ip[[2]]$sumconv[grep('S_msy\\[',conv_f7_ip[[2]]$variable)][f7_ip$median[grep('zstar\\[',f7_ip$variable)]] +
    conv_f7_ip[[2]]$sumconv[grep('zstar\\[',conv_f7_ip[[2]]$variable)]

  Smsy_f8_regime <- apply(drf8[,grep("S_msy\\[",colnames(drf8))],2,postmode)
  Smsy_f8 <- Smsy_f8_regime[zstar_f8]
  Smsy_f8_conv <- 
    conv_f8_ip[[2]]$sumconv[grep('S_msy\\[',conv_f8_ip[[2]]$variable)][f8_ip$median[grep('zstar\\[',f8_ip$variable)]] +
    conv_f8_ip[[2]]$sumconv[grep('zstar\\[',conv_f8_ip[[2]]$variable)]

 

  phmma_Smsy_regime_median <- f6_ip$median[grep('S_msy\\[',f6_ip$variable)]
  phmma_Smsy_median <- phmma_Smsy_regime_median[phmma_zstar_regime_median]
  
  phmmb_Smsy_regime_median<-f7_ip$median[grep('S_msy\\[',f7_ip$variable)]
  phmmb_Smsy_median <-phmmb_Smsy_regime_median[phmmb_zstar_regime_median]

  phmmab_Smsy_regime_median<-f8_ip$median[grep('S_msy\\[',f8_ip$variable)]
  phmmab_Smsy_median <-phmmab_Smsy_regime_median[phmmab_zstar_regime_median]

  smsysim<-smsyCalc(dat$alpha,dat$beta)
  
  dfsmsy<- data.frame(parameter="smsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma","hmmb","hmmab"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(smsysim,8),
                      median=c(rep(f1_ip$median[f1_ip$variable=='S_msy'],nrow(dat)),
                             rep(f2_ip$median[f2_ip$variable=='S_msy'],nrow(dat)),
                             f3_ip$median[grep('S_msy\\[',f3_ip$variable)],
                             f4_ip$median[grep('S_msy\\[',f4_ip$variable)],
                             f5_ip$median[grep('S_msy\\[',f5_ip$variable)],
                             phmma_Smsy_median,
                             phmmb_Smsy_median,
                            phmmab_Smsy_median),
                      mode=c(rep(Smsy_f1,nrow(dat)) ,
                        rep(Smsy_f2,nrow(dat)),
                        Smsy_f3,
                        Smsy_f4,
                        Smsy_f5,
                        Smsy_f6,
                        Smsy_f7,
                        Smsy_f8),
                      convergence=c(rep(conv_f1_ip[[2]]$sumconv[conv_f1_ip[[2]]$variable=="S_msy"],nrow(dat)),
                        rep(conv_f2_ip[[2]]$sumconv[conv_f2_ip[[2]]$variable=="S_msy"],nrow(dat)),
                         conv_f3_ip[[2]]$sumconv[grep('S_msy\\[',conv_f3_ip[[2]]$variable)],
                         conv_f4_ip[[2]]$sumconv[grep('S_msy\\[',conv_f4_ip[[2]]$variable)],
                         conv_f5_ip[[2]]$sumconv[grep('S_msy\\[',conv_f5_ip[[2]]$variable)],
                         Smsy_f6_conv,
                         Smsy_f7_conv,
                         Smsy_f8_conv
                        ),
                        conv_warning=NA )
  
  dfsmsy$pbias<- ((as.numeric(dfsmsy$mode)-dfsmsy$sim)/dfsmsy$sim)*100
  dfsmsy$bias<- (as.numeric(dfsmsy$mode)-dfsmsy$sim)
  
  
  #sgen


  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                       model=rep(c("simple",
                                   "autocorr",
                                   "rwa","rwb","rwab",
                                   "hmma","hmmb","hmmab"),each=nrow(dat)),
                       by=rep(dat$year,8),
                       sim=dat$sGen,
                       median=c(unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="simple"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="simple"], 
                                           b=1/dfsmax$median[dfsmax$model=="simple"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="autocorr"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="autocorr"], 
                                           b=1/dfsmax$median[dfsmax$model=="autocorr"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="rwa"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="rwa"], 
                                           b=1/dfsmax$median[dfsmax$model=="rwa"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="rwb"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="rwb"], 
                                           b=1/dfsmax$median[dfsmax$model=="rwb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="rwab"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="rwab"], 
                                           b=1/dfsmax$median[dfsmax$model=="rwab"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="hmma"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="hmma"], 
                                           b=1/dfsmax$median[dfsmax$model=="hmma"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="hmmb"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="hmmb"], 
                                           b=1/dfsmax$median[dfsmax$model=="hmmb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$median[dfa$model=="hmmab"],
                                           Smsy=dfsmsy$median[dfsmsy$model=="hmmab"], 
                                           b=1/dfsmax$median[dfsmax$model=="hmmab"]))),
                       mode=c(unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="simple"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="simple"], 
                                           b=1/dfsmax$mode[dfsmax$model=="simple"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="autocorr"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="autocorr"], 
                                           b=1/dfsmax$mode[dfsmax$model=="autocorr"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="rwa"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="rwa"], 
                                           b=1/dfsmax$mode[dfsmax$model=="rwa"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="rwb"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="rwb"], 
                                           b=1/dfsmax$mode[dfsmax$model=="rwb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="rwab"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="rwab"], 
                                           b=1/dfsmax$mode[dfsmax$model=="rwab"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="hmma"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="hmma"], 
                                           b=1/dfsmax$mode[dfsmax$model=="hmma"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="hmmb"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="hmmb"], 
                                           b=1/dfsmax$mode[dfsmax$model=="hmmb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$mode[dfa$model=="hmmab"],
                                           Smsy=dfsmsy$mode[dfsmsy$model=="hmmab"], 
                                           b=1/dfsmax$mode[dfsmax$model=="hmmab"]))),
                      convergence=dfsmsy$convergence,
                      conv_warning=NA
                       )
  
  dfsgen$pbias<- ((as.numeric(dfsgen$mode)-dfsgen$sim)/dfsgen$sim)*100   
  dfsgen$bias<- (as.numeric(dfsgen$mode)-dfsgen$sim) 
  
  #umsy
  #S msy
  umsy_f1 <- postmode(x=drf1$"U_msy")
  umsy_f2 <- postmode(x=drf2$"U_msy")

  umsy_f3 <- apply(drf3[,grep("U_msy\\[",colnames(drf3))],2,postmode)
  umsy_f4 <- postmode(x=drf4$"U_msy")
  umsy_f5 <- apply(drf5[,grep("U_msy\\[",colnames(drf5))],2,postmode)
  

  umsy_f6_regime <- apply(drf6[,grep("U_msy\\[",colnames(drf6))],2,postmode)
  umsy_f6 <- umsy_f6_regime[zstar_f6]
  umsy_f6_conv <- 
    conv_f6_ip[[2]]$sumconv[grep('U_msy\\[',conv_f6_ip[[2]]$variable)][f6_ip$median[grep('zstar\\[',f6_ip$variable)]] +
    conv_f6_ip[[2]]$sumconv[grep('zstar\\[',conv_f6_ip[[2]]$variable)]

  umsy_f7 <- postmode(x=drf7$"U_msy")

  umsy_f8_regime <- apply(drf8[,grep("U_msy\\[",colnames(drf8))],2,postmode)
  umsy_f8 <- umsy_f8_regime[zstar_f8]
  umsy_f8_conv <- 
    conv_f8_ip[[2]]$sumconv[grep('U_msy\\[',conv_f8_ip[[2]]$variable)][f8_ip$median[grep('zstar\\[',f8_ip$variable)]] +
    conv_f8_ip[[2]]$sumconv[grep('zstar\\[',conv_f8_ip[[2]]$variable)]

  phmma_umsy_regime_median <- f6_ip$median[grep('U_msy\\[',f6_ip$variable)]
  phmma_umsy_median <- phmma_umsy_regime_median[phmma_zstar_regime_median]
  
  phmmab_umsy_regime_median<-f8_ip$median[grep('U_msy\\[',f8_ip$variable)]
  phmmab_umsy_median <-phmmab_umsy_regime_median[phmmab_zstar_regime_median]


  dfumsy<- data.frame(parameter="umsy",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                      model=rep(c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma","hmmb","hmmab"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(samEst::umsyCalc(dat$alpha),8),
                      median=c(rep(f1_ip$median[f1_ip$variable=='U_msy'],nrow(dat)),
                             rep(f2_ip$median[f2_ip$variable=='U_msy'],nrow(dat)),
                             f3_ip$median[grep('U_msy\\[',f3_ip$variable)],
                             rep(f4_ip$median[f4_ip$variable=='U_msy'],nrow(dat)),
                             f5_ip$median[grep('U_msy\\[',f5_ip$variable)],
                             phmma_umsy_median,
                             rep(f7_ip$median[f7_ip$variable=='U_msy'],nrow(dat)),
                            phmmab_umsy_median),
                      mode=c(rep(umsy_f1,nrow(dat)),
                        rep(umsy_f2,nrow(dat)),
                        umsy_f3,
                        rep(umsy_f4,nrow(dat)),
                        umsy_f5,
                        umsy_f6,
                        rep(umsy_f7,nrow(dat)),
                        umsy_f8),
                      convergence=c(rep(conv_f1_ip[[2]]$sumconv[conv_f1_ip[[2]]$variable=="U_msy"],nrow(dat)),
                        rep(conv_f2_ip[[2]]$sumconv[conv_f2_ip[[2]]$variable=="U_msy"],nrow(dat)),
                         conv_f3_ip[[2]]$sumconv[grep('U_msy\\[',conv_f3_ip[[2]]$variable)],
                         rep(conv_f4_ip[[2]]$sumconv[conv_f4_ip[[2]]$variable=="U_msy"],nrow(dat)),
                         conv_f5_ip[[2]]$sumconv[grep('U_msy\\[',conv_f5_ip[[2]]$variable)],
                         umsy_f6_conv,
                         rep(conv_f7_ip[[2]]$sumconv[conv_f7_ip[[2]]$variable=="U_msy"],nrow(dat)),
                         umsy_f8_conv),
                      conv_warning=NA)
  
  dfumsy$pbias<- ((as.numeric(dfumsy$mode)-dfumsy$sim)/dfumsy$sim)*100
  dfumsy$bias<- (as.numeric(dfumsy$mode)-dfumsy$sim)
  
  ##logliks
  ll=list(f1$draws(variables=c('log_lik'),format='draws_matrix'),
          f2$draws(variables=c('log_lik'),format='draws_matrix'),
          f3$draws(variables=c('log_lik'),format='draws_matrix'),
          f4$draws(variables=c('log_lik'),format='draws_matrix'),
          f5$draws(variables=c('log_lik'),format='draws_matrix'),
          f6$draws(variables=c('log_lik'),format='draws_matrix'),
          f7$draws(variables=c('log_lik'),format='draws_matrix'),
          f8$draws(variables=c('log_lik'),format='draws_matrix'))
  
  #Mean pointwise loglikelihoods
  dfelpd<- data.frame(parameter="mean_ELPD",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",8),
                      model=c("simple",
                                  "autocorr",
                                  "rwa","rwb","rwab",
                                  "hmma","hmmb","hmmab"),
                      by=c(NA,8),
                      sim=rep(NA,8),
                      median=c(sum(apply(ll[[1]],2,log_mean_exp)),
                     sum(apply(ll[[2]],2,samEst::log_mean_exp)),
                    sum(apply(ll[[3]],2,samEst::log_mean_exp)),
                    sum(apply(ll[[4]],2,samEst::log_mean_exp)),
                    sum(apply(ll[[5]],2,samEst::log_mean_exp)),
                    sum(apply(ll[[6]],2,samEst::log_mean_exp)), 
                    sum(apply(ll[[7]],2,samEst::log_mean_exp)),
                    sum(apply(ll[[8]],2,samEst::log_mean_exp))),
                      mode=NA,
                      convergence=rep(NA,8),
                      conv_warning=NA,
                      pbias=rep(NA,8),
                      bias=NA)
  #Stan AIC estimates
  dfaic<- data.frame(parameter="AIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab",
                                 "hmma","hmmb","hmmab"),
                     by=rep(NA,8),
                     sim=rep(NA,8),
                     median=c(samEst::stan_aic(x=ll,form='aic',type='full',k=c(3,4,4,4,5,6,6,7))),
                     mode=NA,
                     convergence=rep(NA,8),
                     conv_warning=NA,
                     pbias=rep(NA,8),
                     bias=NA)
  
  #Stan BIC estimates
  dfbic<- data.frame(parameter="BIC",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                                 "autocorr",
                                 "rwa","rwb","rwab",
                                 "hmma","hmmb","hmmab"),
                     by=rep(NA,8),
                     sim=rep(NA,8),
                     median=c(samEst::stan_aic(x=ll,form='bic',type='full',k=c(3,4,4,4,5,6,6,7))),
                     mode=NA,
                     convergence=rep(NA,8),
                     conv_warning=NA,
                     pbias=rep(NA,8),
                     bias=NA)

  dff<-rbind(dfa,dfsmax,dfsig,dfsiga,dfsigb,dfsmsy,dfsgen,dfumsy,dfelpd,dfaic,dfbic)
  
  return(dff)
  
}




