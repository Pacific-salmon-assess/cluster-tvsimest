

#================================================================================================================
#stan runs
#================================================================================================================

library(cmdstanr)
library(rslurm)
library(samEst)
source("R/stan_rwtest_func.R")


file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f.stan")
mod3=cmdstanr::cmdstan_model(file3)

file3ja=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ja.stan")
mod3_ja=cmdstanr::cmdstan_model(file3ja)

file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)


file4ja=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ja.stan")
mod4_ja=cmdstanr::cmdstan_model(file4ja)


file4jan=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_jan.stan")
mod4_jan=cmdstanr::cmdstan_model(file4jan)


#---------------------------------------------------------------------------------------------------
#stan  base scenarios


simPars <- read.csv("data/generic/SimPars.csv")



pars<-data.frame(path="..",
  a=rep(c(4,5,6,11),each=1000),
  u=1:1000)





sjobstan <- slurm_apply(stan_rwtest_func, pars, jobname = 'stanrwtest',
                    nodes = 120, cpus_per_node = 5, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr"),
                    rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars", "mod3","mod3_ja"
                      "mod4","mod4_ja","mod4_jan"))



#AFTER JOB IS DONE IMPORT  the results
resstan <- get_slurm_out(sjobstan, outtype = 'table', wait = FALSE)

saveRDS(resstan, file = "resstanrw.rds")