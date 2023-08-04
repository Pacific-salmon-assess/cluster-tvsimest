#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# March 2022 
#=============================================================


#remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan")
#
#
# run these if first time running script or if updates were implemented. 

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)

#packageVersion("Matrix")

library(cmdstanr)
library(rslurm)
library(samEst)
source("R/stan_func.R")
source("R/tmb_func.R")
source("R/utils.R")
source("R/check_stan_conv.R")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(rstan)

#source("sgen_functions.R")

#compile stan models
file1=file.path(cmdstanr::cmdstan_path(),'srmodels', "m1f_ip.stan")
mod1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'srmodels', "m2f_ip.stan")
mod2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ip.stan")
mod3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ip.stan")
mod4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f_ip.stan")
mod5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'srmodels', "m6f_ip.stan")
mod6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'srmodels', "m7f_ip.stan")
mod7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'srmodels', "m8f_ip.stan")
mod8=cmdstanr::cmdstan_model(file8)




rstan_options(auto_write = TRUE)
options(mc.cores = 5)


simPars <- read.csv("data/generic/SimPars.csv")


tst_stan<-stan_func(path=".", a=6,u=19)


tst_tmb<-tmb_func(path=".",
  a=6,
  u=19)


head(tst_stan)

head(tst_tmb)

tst<-rbind()

