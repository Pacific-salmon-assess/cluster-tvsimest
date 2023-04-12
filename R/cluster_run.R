library(rslurm)
library(samEst)


simPars <- read.csv("data/generic/SimPars.csv")

source("R/tmb_func.R")

#base case 
tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=5,
  u=1)
  
 
pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 50, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)
head(res, 3)


saveRDS(res, file = "res.rds")

#============================================================================
#sensitivity a scenarios


simPars <- read.csv("data/sensitivity/SimPars.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobtmb_a <- slurm_apply(tmb_func, pars_a, jobname = 'TMBrun',
                    nodes = 50, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))





#---------------------------------------------------------------------------------------------------------
#stan version

library(cmdstanr)

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

source("R/stan_func.R")


tst <- stan_func(path=".",
  a=1,
  u=1)


#I can run a maximum of 6 scenarios at a time, otherwise R cannot run the results in. 

#base case scn 1-6
pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)/2),each=1000),
  u=1:1000)

sjobstan <- slurm_apply(stan_func, pars, jobname = 'stanrun',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8"))



#AFTER JOB IS DONE IMPORT  the results
resstan <- get_slurm_out(sjobstan, outtype = 'table', wait = FALSE)
head(resstan, 3)


saveRDS(resstan, file = "/resstan1_6.rds")

cleanup_files(sjobstan)


#base case scn 7-12
 

pars<-data.frame(path="..",
  a=rep((nrow(simPars)/2):nrow(simPars),each=1000),
  u=1:1000)
 

sjobstan2 <- slurm_apply(stan_func, pars, jobname = 'stanrun',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8"))



#AFTER JOB IS DONE IMPORT  the results
resstan2 <- get_slurm_out(sjobstan2, outtype = 'table', wait = FALSE)
head(resstan2, 3)


saveRDS(resstan, file = "resstan7_12.rds")

cleanup_files(sjobstan)


#log_a sensitivity

#log_a sensitivity - half smax

#Smax sensitivity

#Smax sensitivity

