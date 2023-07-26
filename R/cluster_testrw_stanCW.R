

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
#rw Jacobian adj tests


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



#---------------------------------------------------------------------------------------------------
#informative priors test

library(cmdstanr)
library(rslurm)
library(samEst)
source("R/compare_logbprior_func.R")


file1=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1f.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2f.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f.stan")
mod3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f.stan")
mod5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6f.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7f.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8f.stan")
mod8=cmdstanr::cmdstan_model(file8)


file1_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1f_ip.stan")
mod1_ip=cmdstanr::cmdstan_model(file1_ip)
file2_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2f_ip.stan")
mod2_ip=cmdstanr::cmdstan_model(file2_ip)
file3_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ip.stan")
mod3_ip=cmdstanr::cmdstan_model(file3_ip)
file4_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ip.stan")
mod4_ip=cmdstanr::cmdstan_model(file4_ip)
file5_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f_ip.stan")
mod5_ip=cmdstanr::cmdstan_model(file5_ip)
file6_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6f_ip.stan")
mod6_ip=cmdstanr::cmdstan_model(file6_ip)
file7_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7f_ip.stan")
mod7_ip=cmdstanr::cmdstan_model(file7_ip)
file8_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8f_ip.stan")
mod8_ip=cmdstanr::cmdstan_model(file8_ip)


simPars <- read.csv("data/generic/SimPars.csv")



pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjob_ip <- slurm_apply(compare_logbprior_func, pars, jobname = 'priorcompare',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr"),
                    rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars", "mod3","mod3_ip",
                      "mod4","mod4_ip"))



#AFTER JOB IS DONE IMPORT  the results
res_ip <- get_slurm_out(sjob_ip, outtype = 'table', wait = FALSE)

saveRDS(res_ip, file = "res_ip_sip_compare.rds")


