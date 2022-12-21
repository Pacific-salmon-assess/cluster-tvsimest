library(rslurm)


simPars <- read.csv("data/harcnkSimPars.csv")
save(simPars, file = "data/harcnkSimPars.RData")


test_func <- function(path=".",a, u) {
  
  simData <- list()  
  allsimest <- list()
  simData[[a]] <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  #compiled Bayesian models try moving this out of function
  simple_mod <- samEst::compile_code(type='static', ac=FALSE, par='n')

  dat <- simData[[a]][simData[[a]]$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  p <- ricker_TMB(data=df)

  b <- ricker_stan(data=df,iter = 800, mod=simple_mod)

  dfa<- data.frame(parameter="alpha",
              iteration=u,
              scenario= simPars$scenario[a],
              method=rep(c("MLE","MCMC"),each=nrow(df)),
              model=rep(c("simple",
                "simple"),each=nrow(df)),
              by=rep(dat$year,2),
              sim=rep(dat$alpha,2),
              est=c(rep(p$alpha,nrow(df)),
                    rep(b$alpha,nrow(df))),
              convergence=c(rep(c(p$model$convergence + p$conv_problem,
               as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1)),each=nrow(df))))
  dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
      
  return(dfa)

}

    
pars<-data.frame(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=rep(1,100),
  u=1:100)


sjob <- slurm_apply(test_func, pars, jobname = 'test_apply',
                    nodes = 2, cpus_per_node = 2, submit = FALSE,
                    pkgs=c("samEst","samSim","rstan"),
                    global_objects=c("data/harcnkSimPars.RData"))
