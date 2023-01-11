library(rslurm)
library(samEst,lib="/fs/vnas_Hdfo/comda/caw001/Rlib")
library(rstan,lib="/fs/vnas_Hdfo/comda/caw001/Rlib")


simPars <- read.csv("data/harcnkSimPars.csv")
#save(simPars, file = "data/harcnkSimPars.RData")
#load("data/harcnkSimPars.RData")
simple_mod <- samEst::compile_code(type='static', ac=FALSE, par='n',lambertW = FALSE)
simpleac_mod <- samEst::compile_code(type='static', ac=TRUE, par='n',lambertW = FALSE)
rwa_mod <- samEst::compile_code(type='rw',ac=FALSE,par="a",lambertW = FALSE)
rwb_mod <- samEst::compile_code(type='rw',ac=FALSE,par="b",lambertW = FALSE)
rwab_mod <- samEst::compile_code(type='rw',ac=FALSE,par="both",lambertW = FALSE)
hmma_mod <- samEst::compile_code(type='hmm',ac=FALSE,par="a",lambertW = FALSE)
hmmb_mod <- samEst::compile_code(type='hmm',ac=FALSE,par="b",lambertW = FALSE)
hmmab_mod <- samEst::compile_code(type='hmm',ac=FALSE,par="both",lambertW = FALSE)



test_func <- function(path=".",a, u) {
  
  simData <- list()  
  allsimest <- list()
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


  b <- ricker_stan(data=df,iter = 800, mod=simple_mod)
  bac <- ricker_stan(data=df,iter = 800, AC=TRUE, mod=simpleac_mod)
  btva <- ricker_rw_stan(data=df, par="a",iter = 800, mod=rwa_mod)
  btvb <- ricker_rw_stan(data=df, par="b",iter = 800, mod=rwb_mod)
  btvab <- ricker_rw_stan(data=df, par="both",iter = 800, mod=rwab_mod) 
  bhmma <- ricker_hmm_stan(data=df, par="a",iter = 800, mod=hmma_mod)
  bhmmb <- ricker_hmm_stan(data=df, par="b",iter = 800, mod=hmmb_mod)
  bhmmab <- ricker_hmm_stan(data=df, par="both",iter = 800, mod=hmmab_mod) 
   

  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c(rep("MLE",8),rep("MCMC",8)),each=nrow(df)),
              model=rep(c("simple",
                   "autocorr",
                   "rwa","rwb","rwab",
                   "hmma_regime","hmmb_regime","hmmab_regime",
                   "simple","autocorr",
                   "rwa","rwb","rwab",
                   "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
              by=rep(dat$year,16),
              sim=rep(dat$alpha,16),
              est=c(rep(p$alpha,nrow(df)),
                    rep(pac$alpha,nrow(df)),
                    ptva$alpha,
                    rep(ptvb$alpha,nrow(df)),
                    ptvab$alpha,
                    phmma$alpha[phmma$regime],
                    rep(phmmb$alpha,nrow(df)),
                    phmm$alpha[phmm$regime],
                    rep(b$alpha,nrow(df)),
                    rep(bac$alpha,nrow(df)),
                    btva$alpha[-1],
                    rep(btvb$alpha,nrow(df)),
                    btvab$alpha[-1],
                    bhmma$alpha_regime,
                    rep(bhmmb$alpha,nrow(df)),
                    bhmmab$alpha_regime),
              convergence=c(rep(c(p$model$convergence + p$conv_problem,
                    pac$model$convergence + pac$conv_problem,
                    ptva$model$convergence + ptva$conv_problem,
                    ptvb$model$convergence + ptvb$conv_problem,
                    ptvab$model$convergence + ptvab$conv_problem,
                    phmma$model$convergence + phmma$conv_problem,
                    phmmb$model$convergence + phmmb$conv_problem,
                    phmm$model$convergence + phmm$conv_problem,
                    as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
                    as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1)
                    ),each=nrow(df)),
                    as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
                    rep(as.numeric(abs(btvb$mcmcsummary[grep("log_a",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
                    as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),
                    #hmma pick
                    c(as.numeric(abs(bhmma$mcmcsummary[grep("log_a\\[",rownames(bhmma$mcmcsummary)),
                     "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
                     as.numeric(abs(bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
                    #hmmb 
                    rep(as.numeric(abs(bhmmb$mcmcsummary[grep("log_a",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
                    #hmmab pick
                    c(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a\\[",rownames(bhmmab$mcmcsummary)),
                      "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
                      as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))))
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
      
  return(dfa)

}


tst<-test_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=1,
  u=1)
  
pars<-data.frame(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=rep(1,100),
  u=1:100)


sjob <- slurm_apply(test_func, pars, jobname = 'test_apply',
                    nodes = 1, cpus_per_node = 50, submit = FALSE,
                    pkgs=c("samEst","rstan"),
                    rscript_path = "/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
                    libPaths="/fs/vnas_Hdfo/comda/caw001/Rlib",
                    global_objects=c("simPars",
                      "simple_mod",
                      "simpleac_mod",
                      "rwa_mod", 
                      "rwb_mod", 
                      "rwab_mod", 
                      "hmma_mod", 
                      "hmmb_mod", 
                      "hmmab_mod"))




#modify submit.sh
#which account should I use
#which comment should I use

#https://portal.science.gc.ca/confluence/display/SCIDOCS/Quick+Start+to+Using+Linux+Clusters+With+SLURM