library(rslurm)

library(samEst)


#simPars <- read.csv("data/generic/SimPars.csv")
simPars <- read.csv("data/sensitivity/SimPars.csv")
#save(simPars, file = "data/harcnkSimPars.RData")
#load("data/harcnkSimPars.RData")

source("R/rslurm_TMB_func.R")



tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=1,
  u=1)
  R
pars<-data.frame(path="..",
  a=rep(1:8,each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 1, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
                    libPaths="/fs/vnas_Hdfo/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res <- get_slurm_out(sjobtmb, outtype = 'table', wait = FALSE)
head(res, 3)


#https://portal.science.gc.ca/confluence/display/SCIDOCS/Quick+Start+to+Using+Linux+Clusters+With+SLURM

cmdstanpth <- read

#Stan
cmdstanr::set_cmdstan_path(path="/fs/vnas_Hdfo/comda/dag004/.cmdstan/cmdstan-2.31.0")
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1f.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f.stan")
mod3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5f.stan")
mod5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6f.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7f.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8f.stan")
mod8=cmdstanr::cmdstan_model(file8)
