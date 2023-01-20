library(rslurm)
library(samEst,lib="/fs/vnas_Hdfo/comda/caw001/Rlib")


simPars <- read.csv("data/harcnkSimPars.csv")
#save(simPars, file = "data/harcnkSimPars.RData")
#load("data/harcnkSimPars.RData")



tmb_func <- function(path=".",a, u) {
  
  simData <- list()  

  simData[[a]] <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  #compiled Bayesian models try moving this out of function
  
  dat <- simData[[a]][simData[[a]]$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))


  p <- ricker_TMB(data=df)
  pac <- ricker_TMB(data=df, AC=TRUE)
  ptva <- ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- ricker_rw_TMB(data=df, tv.par='both')
  phmma <- ricker_hmm_TMB(data=df, tv.par='a')
  phmmb <- ricker_hmm_TMB(data=df, tv.par='b')
  phmm  <- ricker_hmm_TMB(data=df, tv.par='both')



  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",8)),each=nrow(df)),
              model=rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
              by=rep(dat$year,8),
              sim=rep(dat$alpha,8),
              est=c(rep(p$alpha,nrow(df)),
                    rep(pac$alpha,nrow(df)),
                    ptva$alpha,
                    rep(ptvb$alpha,nrow(df)),
                    ptvab$alpha,
                    phmma$alpha[phmma$regime],
                    rep(phmmb$alpha,nrow(df)),
                    phmm$alpha[phmm$regime]),
              convergence=c(rep(c(p$model$convergence + p$conv_problem,
                    pac$model$convergence + pac$conv_problem,
                    ptva$model$convergence + ptva$conv_problem,
                    ptvb$model$convergence + ptvb$conv_problem,
                    ptvab$model$convergence + ptvab$conv_problem,
                    phmma$model$convergence + phmma$conv_problem,
                    phmmb$model$convergence + phmmb$conv_problem,
                    phmm$model$convergence + phmm$conv_problem
                    ),each=nrow(df))))
                    
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
      
  return(dfa)

}


tst<-tmb_func(path=".",
  a=1,
  u=1)
  
pars<-data.frame(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=rep(1:12,each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 1, cpus_per_node = 50, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
                    libPaths="/fs/vnas_Hdfo/comda/caw001/Rlib",
                    global_objects=c("simPars"))




res <- get_slurm_out(sjobtmb, outtype = 'table', wait = FALSE)
head(res, 3)


#https://portal.science.gc.ca/confluence/display/SCIDOCS/Quick+Start+to+Using+Linux+Clusters+With+SLURM