
##remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)
#================================================================================================================
#stan runs
#================================================================================================================

library(cmdstanr)
library(rslurm)
library(samEst)
source("R/stan_func.R")
source("R/utils.R")
source("R/check_stan_conv.R")


file1=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1f_ip.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2f_ip.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ip.stan")
mod3=cmdstanr::cmdstan_model(file3)

file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_smax.stan")
mod4=cmdstanr::cmdstan_model(file4)
#file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ip.stan")

file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f_smax.stan")
mod5=cmdstanr::cmdstan_model(file5)
#file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f_ip.stan")


file6=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6f_ip.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7f_ip.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8f_ip.stan")
mod8=cmdstanr::cmdstan_model(file8)



#---------------------------------------------------------------------------------------------------
#stan  base scenarios


simPars <- read.csv("data/generic/SimPars.csv")
# a=1
# u=40

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=200),
  u=1:200)

#run on Aug 1
#pars<-data.frame(path="..",
#  a=rep(c(1,2),each=1000),
#  u=1:1000)

#pars<-data.frame(path="..",
#  a=rep(seq_len(nrow(simPars)),each=10),
#  u=1:10)

tst<-stan_func(path=".", a=6,u=19)

sjobstan <- slurm_apply(stan_func, pars, jobname = 'stanrunhi',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr"),
                    rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode", "check_stan_conv" ))



#AFTER JOB IS DONE IMPORT  the results
resstan <- get_slurm_out(sjobstan, outtype = 'table', wait = FALSE)

saveRDS(resstan, file = "resstan.rds")
saveRDS(resstan[resstan$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan1.rds")
saveRDS(resstan[resstan$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan2.rds")

#---------------------------------------------------------------------------------------------------
#stan  sensitivity a scenarios
#(export TMPDIR="/home/caw001/Documents/tvsimest/stantmp" ; R )


simPars <- read.csv("data/sensitivity/SimPars.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstan_a <- slurm_apply(stan_func, pars_a, jobname = 'stanrun_a',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))




#AFTER JOB IS DONE IMPORT  the results
resstan_a <- get_slurm_out(sjobstan_a, outtype = 'table', wait = TRUE)

head(resstan_a, 3)

saveRDS(resstan_a[resstan_a$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan_a1.rds")
saveRDS(resstan_a[resstan_a$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan_a2.rds")
saveRDS(resstan_a, file = "resstan_a.rds")


#---------------------------------------------------------------------------------------------------
#stan sensitivity a halfsmax scenarios

simPars <- read.csv("data/sensitivity_halfSmax/SimPars.csv")

pars_asmax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstan_asmax <- slurm_apply(stan_func, pars_asmax, jobname = 'stanrun_asmax',
                            nodes = 250, cpus_per_node = 5, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))




#AFTER JOB IS DONE IMPORT  the results
resstan_asmax <- get_slurm_out(sjobstan_asmax, outtype = 'table', wait = TRUE)

head(resstan_asmax, 3)

saveRDS(resstan_asmax[resstan_asmax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan_asmax1.rds")
saveRDS(resstan_asmax[resstan_asmax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan_asmax2.rds")
saveRDS(resstan_asmax, file = "resstan_asmax.rds")

#---------------------------------------------------------------------------------------------------
#smax scenarios 




simPars <- read.csv("data/Smax_sensitivity/SimPars.csv")

pars_smax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstan_smax <- slurm_apply(stan_func, pars_smax, jobname = 'stanrun_smax',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))




#AFTER JOB IS DONE IMPORT  the results
resstan_smax <- get_slurm_out(sjobstan_smax, outtype = 'table', wait = TRUE)

head(resstan_smax, 3)

saveRDS(resstan_smax[resstan_smax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan_smax1.rds")
saveRDS(resstan_smax[resstan_smax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan_smax2.rds")
saveRDS(resstan_smax, file = "resstan_smax.rds")


#---------------------------------------------------------------------------------------------------
#smax scenarios double alpha

simPars <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")

pars_smaxda<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstan_smaxda <- slurm_apply(stan_func, pars_smaxda, jobname = 'stanrun_smaxda',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8"))




#AFTER JOB IS DONE IMPORT  the results
resstan_smaxda <- get_slurm_out(sjobstan_smaxda, outtype = 'table', wait = TRUE)

head(resstan_smaxda, 3)

saveRDS(resstan_smaxda[resstan_smaxda$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan_smaxda1.rds")
saveRDS(resstan_smaxda[resstan_smaxda$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan_smaxda2.rds")
saveRDS(resstan_smaxda, file = "resstan_smaxda.rds")




#---------------------------------------------------------------------------------------------------
#sensitivity sigma scenarios low

simPars <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

pars_siglow<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstan_siglow <- slurm_apply(stan_func, pars_siglow, jobname = 'stanrun_siglow',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))


resstan_siglow <- get_slurm_out(sjobstan_siglow, outtype = 'table', wait = TRUE)



head(resstan_siglow, 3)

saveRDS(resstan_siglow[resstan_siglow$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan_siglow1.rds")
saveRDS(resstan_siglow[resstan_siglow$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan_siglow2.rds")
saveRDS(resstan_siglow, file = "resstan_siglow.rds")





#---------------------------------------------------------------------------------------------------
#sensitivity sigma scenarios med

simPars <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

pars_sigmed<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstan_sigmed <- slurm_apply(stan_func, pars_sigmed, jobname = 'stanrun_sigmed',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))


resstan_sigmed <- get_slurm_out(sjobstan_sigmed, outtype = 'table', wait = TRUE)


saveRDS(resstan_sigmed[resstan_sigmed$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "resstan_sigmed1.rds")
saveRDS(resstan_sigmed[resstan_sigmed$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "resstan_sigmed2.rds")
saveRDS(resstan_sigmed, file = "resstan_sigmed.rds")




#---------------------------------------------------------------------------------------------------
#Base ER

simPars <- read.csv("data/genericER/SimPars_ER.csv")

pars_baseER<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstan_baseER<- slurm_apply(stan_func, pars_baseER, jobname = 'stanrun_baseER',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8","postmode","check_stan_conv"))


resstan_baseER <- get_slurm_out(sjobstan_baseER, outtype = 'table', wait = TRUE)


saveRDS(resstan_baseER[resstan_baseER$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "resstan_baseER1.rds")
saveRDS(resstan_baseER[resstan_baseER$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "resstan_baseER2.rds")
saveRDS(resstan_baseER, file = "resstan_baseER.rds")



#---------------------------------------------------------------------------------------------------
#No need to run bias corr for stan as it does not make much sense in a bayesian context.





#================================================================================================================
#loo runs
#================================================================================================================


#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)
#base
library(rslurm)
library(samEst)
source("R/stan_lfo_func.R") #stan lfo function
library(cmdstanr)



#load in cmdstanr models for LFO
file1=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1loo_ip.stan")
mod1lfo=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2loo_ip.stan")
mod2lfo=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3loo_ip.stan")
mod3lfo=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4loo_smax_ip.stan")
mod4lfo=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5loo_smax_ip.stan")
mod5lfo=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6loo_ip.stan")
mod6lfo=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7loo_ip.stan")
mod7lfo=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8loo_ip.stan")
mod8lfo=cmdstanr::cmdstan_model(file8)



#simulation parameters
simPars <- read.csv("data/generic/SimPars.csv")

#test function- this takes awhile....
#tst <- stan_lfo(path=".",
#                a=1,
#                u=1)

#test some pars

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



#pars<-data.frame(path="..",
#  a=rep(2,each=15),
#  u=1:15)
#
#stan_lfo(path=".", a=4,u=23)

#slurm job
sjobstanloobase <- slurm_apply(stan_lfo, pars, jobname = 'stanloobase',
                            nodes = 300, cpus_per_node = 3, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo <- get_slurm_out(sjobstanloobase, outtype = 'table', wait = FALSE)

saveRDS(resstanloo, file = "resstanloo.rds")
saveRDS(resstanloo[resstanloo$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo1.rds")
saveRDS(resstanloo[resstanloo$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo2.rds")


#---------------------------------------------------------------------------------------------------
#stan loo sensitivity a scenarios
library(rslurm)
library(samEst)
source("R/stan_lfo_func.R") #stan lfo function
library(cmdstanr)


simPars <- read.csv("data/sensitivity/SimPars.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstanloo_a <- slurm_apply(stan_lfo, pars_a, jobname = 'stanloo_a',
                            nodes = 250, cpus_per_node = 3, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))




#AFTER JOB IS DONE IMPORT  the results
resstanloo_a <- get_slurm_out(sjobstanloo_a, outtype = 'table', wait = TRUE)

head(resstan_a, 3)

saveRDS(resstanloo_a[resstanloo_a$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_a1.rds")
saveRDS(resstanloo_a[resstanloo_a$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_a2.rds")
saveRDS(resstanloo_a, file = "resstanloo_a.rds")


#---------------------------------------------------------------------------------------------------
#sensitivity a scenarios half smax
library(rslurm)
library(samEst)
source("R/stan_lfo_func.R") #stan lfo function
library(cmdstanr)

simPars <- read.csv("data/sensitivity_halfSmax/SimPars.csv")

pars_asmax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstanloo_asmax <- slurm_apply(stan_lfo, pars_asmax, jobname = 'stanloo_asmax',
                            nodes = 250, cpus_per_node = 3, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_asmax <- get_slurm_out(sjobstanloo_asmax, outtype = 'table', wait = TRUE)

head(resstanloo_asmax, 3)

saveRDS(resstanloo_asmax[resstanloo_asmax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_asmax1.rds")
saveRDS(resstanloo_asmax[resstanloo_asmax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_asmax2.rds")
saveRDS(resstanloo_asmax, file = "resstanloo_asmax.rds")

#---------------------------------------------------------------------------------------------------
#smax scenarios 
library(rslurm)
library(samEst)
library(cmdstanr)
source("R/stan_lfo_func.R") #stan lfo function



simPars <- read.csv("data/Smax_sensitivity/SimPars.csv")

pars_smax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstanloo_smax <- slurm_apply(stan_lfo, pars_smax, jobname = 'stanloo_smax',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_smax <- get_slurm_out(sjobstanloo_smax, outtype = 'table', wait = TRUE)

head(resstanloo_smax, 3)

saveRDS(resstanloo_smax[resstanloo_smax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_smax1.rds")
saveRDS(resstanloo_smax[resstanloo_smax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_smax2.rds")
saveRDS(resstanloo_smax, file = "resstanloo_smax.rds")


#---------------------------------------------------------------------------------------------------
#smax scenarios double alpha
library(rslurm)
library(samEst)
library(cmdstanr)
source("R/stan_lfo_func.R") #stan lfo function


simPars <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")

pars_smaxda<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)




sjobstanloo_smaxda <- slurm_apply(stan_lfo, pars_smaxda, jobname = 'stanloo_smaxda',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_smaxda <- get_slurm_out(sjobstanloo_smaxda, outtype = 'table', wait = TRUE)



head(resstanloo_smaxda, 3)

saveRDS(resstanloo_smaxda[resstanloo_smaxda$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_smaxda1.rds")
saveRDS(resstanloo_smaxda[resstanloo_smaxda$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_smaxda2.rds")
saveRDS(resstanloo_smaxda, file = "resstanloo_smaxda.rds")



#---------------------------------------------------------------------------------------------------
#sensitivity sigma scenarios low
library(rslurm)
library(samEst)
library(cmdstanr)
source("R/stan_lfo_func.R") #stan lfo function


simPars <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

pars_siglow<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)



sjobstanloo_siglow <- slurm_apply(stan_lfo, pars_siglow, jobname = 'stanloo_siglow',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_siglow <- get_slurm_out(sjobstanloo_siglow, outtype = 'table', wait = TRUE)



head(resstanloo_siglow, 3)

saveRDS(resstanloo_siglow[resstanloo_siglow$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_siglow1.rds")
saveRDS(resstanloo_siglow[resstanloo_siglow$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_siglow2.rds")
saveRDS(resstanloo_siglow, file = "resstanloo_siglow.rds")

#---------------------------------------------------------------------------------------------------
#sensitivity sigma scenarios med
library(rslurm)
library(samEst)
library(cmdstanr)
source("R/stan_lfo_func.R") #stan lfo function

simPars <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

pars_sigmed<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)

sjobstanloo_sigmed <- slurm_apply(stan_lfo, pars_sigmed, jobname = 'stanloo_sigmed',
                            nodes = 250, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_sigmed <- get_slurm_out(sjobstanloo_sigmed, outtype = 'table', wait = TRUE)

head(resstanloo_sigmed, 3)

saveRDS(resstanloo_sigmed[resstanloo_sigmed$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstanloo_sigmed1.rds")
saveRDS(resstanloo_sigmed[resstanloo_sigmed$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstanloo_sigmed2.rds")
saveRDS(resstanloo_sigmed, file = "resstanloo_sigmed.rds")





#---------------------------------------------------------------------------------------------------
#Base ER

simPars <- read.csv("data/genericER/SimPars_ER.csv")

pars_baseER<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobstanloo_baseER <- slurm_apply(stan_lfo, pars_baseER, jobname = 'stanloo_baseER',
                            nodes = 250, cpus_per_node = 3, submit = FALSE,
                            pkgs=c("cmdstanr", "samEst"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/caw001/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))


resstanloo_baseER <- get_slurm_out(sjobstanloo_baseER, outtype = 'table', wait = TRUE)

head(resstanloo_baseER, 3)

saveRDS(resstanloo_baseER[resstanloo_baseER$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "resstanloo_baseER1.rds")
saveRDS(resstanloo_baseER[resstanloo_baseER$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "resstanloo_baseER2.rds")
saveRDS(resstanloo_baseER, file = "resstanloo_baseER.rds")



