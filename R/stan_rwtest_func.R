stan_rwtest_func<- function(path=".", a,u){
  
  
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
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2)
  )


  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  
  #create folder to hold temp files
#  dir.create(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u,sep=''))
  
#  ls=list.files("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/")
#  if(length(ls)>5){
 #   unlink(paste("/fs/vnas_Hdfo/comda/dag004/Rlib/tmp/tmp_cmdst/",u-2,
#                 '/*',sep=''))
#}
  
  options(mc.cores = 5)
  #
  
  
  f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf3noadj<-f3$summary()

  f3_ja <- mod3_ja$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf3ja<-f3_ja$summary()
  
  f4 <- mod4$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4noadj<-f4$summary()
  

  f4_ja <- mod4_ja$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4ja<-f4_ja$summary()
 

  f4_jan <- mod4_jan$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.95,
                    max_treedepth = 15)
  resf4jan<-f4_jan$summary()

  ptva <- ricker_rw_TMB(data=df_tmb,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df_tmb, tv.par='b')
  #Max. prod
  
  #regime state sequence:
  
  
  dfa<- data.frame(parameter="alpha",
                   iteration=u,
                   scenario= simPars$scenario[a],
                   method=rep(c(rep("MCMC",5),rep("MLE",2)),each=nrow(dat)),
                   model=rep(c("rwa","rwa_ja","rwb","rwb_ja","rwb_jan","rwa_tmb","rwb_tmb"),each=nrow(dat)),
                   by=rep(dat$year,7),
                   sim=rep(dat$alpha,7),
                   est=c(resf3noadj[grep("log_a\\[",resf3noadj$variable),"median"][[1]],
                        resf3ja[grep("log_a\\[",resf3ja$variable),"median"][[1]],
                        rep(resf4noadj[grep("log_a",resf4noadj$variable),"median"][[1]],nrow(dat)),
                        rep(resf4ja[grep("log_a",resf4ja$variable),"median"][[1]],nrow(dat)),
                        rep(resf4jan[grep("log_a",resf4jan$variable),"median"][[1]],nrow(dat)),
                        ptva$alpha,
                        rep(ptvb$alpha,nrow(dat))
                         ),
                   convergence= as.numeric( c(
                        abs(resf3noadj[grep("log_a\\[",resf3noadj$variable),"rhat"][[1]]-1)>.1,
                        abs(resf3ja[grep("log_a\\[",resf3ja$variable),"rhat"][[1]]-1)>.1,
                        rep(abs(resf4noadj[grep("log_a",resf4noadj$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                        rep(abs(resf4ja[grep("log_a",resf4ja$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                        rep(abs(resf4jan[grep("log_a",resf4jan$variable),"rhat"][[1]]-1)>.1,nrow(dat)),

                        rep(ptva$model$convergence + ptva$conv_problem,nrow(dat)),
                        rep(ptvb$model$convergence + ptvb$conv_problem,nrow(dat)))))


  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  dfa$bias<- (dfa$est-dfa$sim)
 
  #Smax
 
  dfsmax<- data.frame(parameter="smax",
                      iteration=u,
                      scenario= simPars$scenario[a],
                      method=rep(c(rep("MCMC",5),rep("MCMC",2)),each=nrow(dat)),
                      model=rep(c("rwa","rwa_ja","rwb","rwb_ja","rwb_jan","rwa_tmb","rwb_tmb"),each=nrow(dat)),
                      by=rep(dat$year,7),
                      sim=rep(1/dat$beta,7),
                      est=c(rep(resf3noadj[grep('S_max',resf3noadj$variable),"median"][[1]],nrow(dat)),
                            rep(resf3ja[grep('S_max',resf3ja$variable),"median"][[1]],nrow(dat)),
                            resf4noadj[grep('S_max',resf4noadj$variable),"median"][[1]],
                            resf4ja[grep('S_max',resf4ja$variable),"median"][[1]],
                            resf4jan[grep('S_max',resf4jan$variable),"median"][[1]],
                            rep(ptva$Smax,nrow(dat)),
                            ptvb$Smax),
                      convergence=c(rep(abs(resf3noadj[grep('S_max',resf3noadj$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                            rep(abs(resf3ja[grep('S_max',resf3ja$variable),"rhat"][[1]]-1)>.1,nrow(dat)),
                            abs(resf4noadj[grep('S_max',resf4noadj$variable),"rhat"][[1]]-1)>.1,
                            abs(resf4ja[grep('S_max',resf4ja$variable),"rhat"][[1]]-1)>.1,
                            abs(resf4jan[grep('S_max',resf4jan$variable),"rhat"][[1]]-1)>.1,
                            rep(ptva$model$convergence + ptva$conv_problem,nrow(dat)),
                            rep(ptvb$model$convergence + ptvb$conv_problem,nrow(dat))
                            ))
            
  dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100
  dfsmax$bias<- (dfsmax$est-dfsmax$sim)
  
  
  dff<-rbind(dfa,dfsmax)
  
  return(dff)
  
}
