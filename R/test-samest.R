


install.packages("TMB",lib="/fs/vnas_Hdfo/comda/caw001/Rlib")
install.packages("gsl",lib="/fs/vnas_Hdfo/comda/caw001/Rlib")

#install.packages("rstan")
remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan", lib="/fs/vnas_Hdfo/comda/caw001/Rlib")
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE, lib="/fs/vnas_Hdfo/comda/caw001/Rlib")


library(TMB)
library(rstan)
library(samEst)


df <- read.csv("data/examplesr.csv")


compile("src/Ricker_simple.cpp")
dyn.load(dynlib("src/Ricker_simple"))



parameters_simple<- list(
    alpha=1.5,
    logbeta = log(1e-08),
    logsigobs=log(.4)
    )

SRdata<-list(obs_S=df$S,
    obs_logRS=df$logRS,
    priors=1)

obj_simple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")
  newtonOption(obj_simple, smartsearch=FALSE)

  opt_simple <- nlminb(obj_simple$par,obj_simple$fn,obj_simple$gr)
  rep_simple <- obj_simple$report()



 datm = list(N=nrow(df),
                R_S =df$logRS,
                S=df$S)

fit1 <- stan(
  file = "src/ricker_linear.stan",  # Stan program
  data = datm,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)

  )


dat <- data.frame(by=df$by,
                  S=df$S,
                  R=df$r,
                  logRS=df$logRS)


p <- ricker_TMB(data=dat)

b <- ricker_stan(data=dat,iter = 800, mod=simple_mod)
