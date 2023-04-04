


stan_func<- function(path=".", a,u){
  
  
  allsimest <- list()
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
  phmma_alpha_conv= abs(f6$summary(variables=c('log_a'))$rhat-1)>.1
  phmma_alpha_conv_regime=as.numeric(phmma_alpha_conv[f6$summary(variables=c('zstar'),'median')$median])

  phmmab_alpha=f8$summary(variables=c('log_a'),'median')$median
  phmmab_alpha_regime=phmmab_alpha[f8$summary(variables=c('zstar'),'median')$median]
  phmmab_alpha_conv=abs(f8$summary(variables=c('log_a'))$rhat-1)>.1
  phmmab_alpha_conv_regime=as.numeric(phmmab_alpha_conv[f8$summary(variables=c('zstar'),'median')$median])
  
  


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
                   est=c(rep(f1$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         rep(f2$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f3$summary(variables=c('log_a'),'median')$median,
                         rep(f4$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         f5$summary(variables=c('log_a'),'median')$median,
                         phmma_alpha_regime,
                         rep(f7$summary(variables=c('log_a'),'median')$median,nrow(dat)),
                         phmmab_alpha_regime),
                   convergence= as.numeric(c(rep(abs(f1$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           rep(abs(f2$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           abs(f3$summary(variables=c('log_a'))$rhat-1)>.1,
                           rep(abs(f4$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           abs(f5$summary(variables=c('log_a'))$rhat-1)>.1,
                           phmma_alpha_conv_regime,
                           rep(abs(f7$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           phmmab_alpha_conv_regime)))

  
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
 
  #Smax
  phmmb_smax=f7$summary(variables=c('S_max'),'median')$median
  phmmb_smax_regime=phmmb_smax[f7$summary(variables=c('zstar'),'median')$median]
  phmmb_smax_conv=abs(f7$summary(variables=c('S_max'))$rhat-1)>.1
  phmmb_smax_conv_regime=as.numeric(phmmb_smax_conv[f7$summary(variables=c('zstar'),'median')$median])


  phmmab_smax=f8$summary(variables=c('S_max'),'median')$median
  phmmab_smax_regime=phmmb_smax[f8$summary(variables=c('zstar'),'median')$median]
  phmmab_smax_conv=abs(f8$summary(variables=c('S_max'))$rhat-1)>.1
  phmmab_smax_conv_regime=as.numeric(phmmab_smax_conv[f8$summary(variables=c('zstar'),'median')$median])
  
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
                      est=c(rep(f1$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f2$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            rep(f3$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            f4$summary(variables=c('S_max'),'median')$median,
                            f5$summary(variables=c('S_max'),'median')$median,
                            rep(f6$summary(variables=c('S_max'),'median')$median,nrow(dat)),
                            phmmb_smax_regime,
                            phmmab_smax_regime),
                      convergence=as.numeric(c(rep(abs(f1$summary(variables=c('S_max'))$rhat-1)>.1,nrow(dat)),
                           rep(abs(f2$summary(variables=c('S_max'))$rhat-1)>.1,nrow(dat)),
                           rep(abs(f3$summary(variables=c('S_max'))$rhat-1)>.1,nrow(dat)),
                           rep(abs(f4$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           abs(f5$summary(variables=c('log_a'))$rhat-1)>.1,
                           phmma_alpha_conv_regime,
                           rep(abs(f7$summary(variables=c('log_a'))$rhat-1)>.1,nrow(dat)),
                           phmmab_alpha_conv_regime)))
                      
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  
  #obs error
  dfsig<- data.frame(parameter="sigma",
                     iteration=u,
                     scenario= simPars$scenario[a],
                     method=rep("MCMC",8),
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "hmma","hmmb","hmmab"),
                     by=NA,
                     sim=rep(dat$sigma,8),
                     est=c(f1$summary(variables=c('sigma'),'median')$median,
                           f2$summary(variables=c('sigma'),'median')$median,
                           f3$summary(variables=c('sigma'),'median')$median,
                           f4$summary(variables=c('sigma'),'median')$median,
                           f5$summary(variables=c('sigma'),'median')$median,
                           f6$summary(variables=c('sigma'),'median')$median,
                           f7$summary(variables=c('sigma'),'median')$median,
                           f8$summary(variables=c('sigma'),'median')$median),
                    convergence= as.numeric(c(abs(f1$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f2$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f3$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f4$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f5$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f6$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f7$summary(variables=c('sigma'))$rhat-1)>.1,
                      abs(f8$summary(variables=c('sigma'))$rhat-1)>.1)))
  
  dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
  #sigma a
  dfsiga<- data.frame(parameter="sigma_a",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep("MCMC",2),
                      model=c("rwa","rwab"),
                      by=NA,
                      sim=NA,
                      est=c(f3$summary(variables=c('sigma_a'),'median')$median,
                            f5$summary(variables=c('sigma_a'),'median')$median),
                      convergence=as.numeric(c(abs(f3$summary(variables=c('sigma_a'))$rhat-1)>.1,
                            abs(f5$summary(variables=c('sigma_a'))$rhat-1)>.1)))
  
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
                            f5$summary(variables=c('sigma_b'),'median')$median),
                      convergence=as.numeric(c(abs(f4$summary(variables=c('sigma_b'))$rhat-1)>.1,
                            abs(f5$summary(variables=c('sigma_b'))$rhat-1)>.1)))
  
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
                                  "hmma","hmmb","hmmab"),each=nrow(dat)),
                      by=rep(dat$year,8),
                      sim=rep(smsysim,8),
                      est=c(samEst::smsyCalc(a=dfa$est[dfa$model=="simple"],b=1/dfsmax$est[dfsmax$model=="simple"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="autocorr"],b=1/dfsmax$est[dfsmax$model=="autocorr"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwa"],b=1/dfsmax$est[dfsmax$model=="rwa"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwb"],b=1/dfsmax$est[dfsmax$model=="rwb"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="rwab"],b=1/dfsmax$est[dfsmax$model=="rwab"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmma"],b=1/dfsmax$est[dfsmax$model=="hmma"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmb"],b=1/dfsmax$est[dfsmax$model=="hmmb"]),
                            samEst::smsyCalc(a=dfa$est[dfa$model=="hmmab"],b=1/dfsmax$est[dfsmax$model=="hmmab"])),
                      convergence=dfa$convergence+dfsmax$convergence
                      )
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100
  
  
  dfsgen <- data.frame(parameter="sgen",
                       iteration=u,
                       scenario= simPars$scenario[a],
                       method=rep(c(rep("MCMC",8)),each=nrow(dat)),
                       model=rep(c("simple",
                                   "autocorr",
                                   "rwa","rwb","rwab",
                                   "hmma","hmmb","hmmab"),each=nrow(dat)),
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
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmma"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmma"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmma"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmb"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmb"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmb"])),
                             unlist(mapply(samEst::sGenCalc,a=dfa$est[dfa$model=="hmmab"],
                                           Smsy=dfsmsy$est[dfsmsy$model=="hmmab"], 
                                           b=1/dfsmax$est[dfsmax$model=="hmmab"]))),
                      convergence=dfa$convergence+dfsmax$convergence
                       )
  
  dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100     
  
  #umsy
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
                      est=c(samEst::umsyCalc(dfa$est[dfa$model=="simple"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="autocorr"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwa"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwb"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="rwab"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmma"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmb"]),
                            samEst::umsyCalc(dfa$est[dfa$model=="hmmab"])),
                      convergence=dfa$convergence
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
                                  "hmma","hmmb","hmmab"),each=nrow(dat)),
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
                        f8$summary(variables=c('log_lik'),'median')$median),
                      convergence=rep(NA,8),
                      pbias=rep(NA,8))
  
  
  dff<-rbind(dfa,dfsmax,dfsig,dfsiga,dfsigb,dfsmsy,dfsgen,dfumsy,dfelpd)
  
  return(dff)
  
}
