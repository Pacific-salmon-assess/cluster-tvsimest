

#============================================================================
#sbase case scenarios

library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/generic/SimPars.csv")


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



res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)
#res <- my_get_slurm_out(slr_job_name='TMBrun', nodes.list=0:39, outtype = 'table')





#AFTER JOB IS DONE IMPORT  the results

saveRDS(res[res$scenario%in%simPars$scenario[1:6],], file = "res1.rds")
saveRDS(res[res$scenario%in%simPars$scenario[7:12],], file = "res2.rds")
saveRDS(res, file = "res.rds")


tmb_func(path=".",
  a=5,
  u=1)
#run 14 that failed
.rslurm_id <- 14
.rslurm_istart <- (.rslurm_id)* 240 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 240, nrow(pars))
rslurm_result14<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_result14[[i-3360]]<-tmb_func(path=".",
  a=pars$a[i],
  u=pars$u[i])
}
result14<-do.call(rbind, rslurm_result14)
saveRDS(result14, file = "res14.rds")


#============================================================================
#sensitivity a scenarios
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sensitivity/SimPars.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_a <- slurm_apply(tmb_func, pars_a, jobname = 'TMBrun_a',
                    nodes = 200, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res_a <- get_slurm_out(sjobtmb_a, outtype = 'table', wait = TRUE)

head(res_a, 3)

saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_a1.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_a2.rds")
saveRDS(res_a, file = "res_a.rds")



tmb_func(path=".",
  a=5,
  u=1)
#run 58 that failed
.rslurm_id <- 1
.rslurm_istart <- (.rslurm_id)* 50 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 50, nrow(pars_asmax))
rslurm_result1<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
  print(i)
   rslurm_result1[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_a$a[i],
  u=pars_a$u[i])
}
result1<-do.call(rbind, rslurm_result1)
saveRDS(result1, file = "res_a1.rds")

#The following files are missing: results_1.RDS


#============================================================================
#sensitivity a scenarios half smax
library(rslurm)
library(samEst)
source("R/tmb_func.R")

simPars <- read.csv("data/sensitivity_halfSmax/SimPars.csv")

pars_asmax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_asmax <- slurm_apply(tmb_func, pars_asmax, jobname = 'TMBrun_asmax',
                    nodes = 100, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))

res_asmax <- get_slurm_out(sjobtmb_asmax, outtype = 'table', wait = TRUE)

head(res_asmax, 3)


saveRDS(res_asmax[res_asmax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_asmax1.rds")
saveRDS(res_asmax[res_asmax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_asmax2.rds")
saveRDS(res_asmax, file = "res_asmax.rds")



tmb_func(path=".",
  a=5,
  u=1)
#run 58 that failed
.rslurm_id <- 58
.rslurm_istart <- (.rslurm_id)* 100 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 100, nrow(pars_asmax))
rslurm_result58<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_result58[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_asmax$a[i],
  u=pars_asmax$u[i])
}
result58<-do.call(rbind, rslurm_result58)
saveRDS(result58, file = "res_asmax58.rds")




#============================================================================
#smax scenarios 
library(rslurm)
library(samEst)
source("R/tmb_func.R")

simPars <- read.csv("data/Smax_sensitivity/SimPars.csv")

pars_smax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_smax <- slurm_apply(tmb_func, pars_smax, jobname = 'TMBrun_smax',
                    nodes = 100, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))


res_smax <- get_slurm_out(sjobtmb_smax, outtype = 'table', wait = TRUE)


res_smax <- my_get_slurm_out('TMBrun_smax', nodes.list=0:99, outtype = "table") 
head(res_smax, 3)

saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_smax1.rds")
saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_smax2.rds")
saveRDS(res_smax, file = "res_smax.rds")

#run 95 that failed
.rslurm_id <- 95
.rslurm_istart <- (.rslurm_id)* 100 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 100, nrow(pars_smax))
rslurm_res_smax95<-list()q()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_smax95[[i-9500]]<-tmb_func(path=".",
  a=pars_smax$a[i],
  u=pars_smax$u[i])
}
result_smax_95<-do.call(rbind, rslurm_res_smax95)
saveRDS(result_smax_95, file = "res_smax_95.rds")


#============================================================================
#smax scenarios double alpha
library(rslurm)
library(samEst)
source("R/tmb_func.R")

simPars <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")

pars_smaxda<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)

#test
paste0(".","/outs/SamSimOutputs/simData/", simPars$nameOM[5],"/",simPars$scenario[5],"/",
                         paste(simPars$nameOM[5],"_", simPars$nameMP[5], "_", "CUsrDat.RData",sep=""))
tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
  a=5,
  u=1)
  

sjobtmb_smaxda <- slurm_apply(tmb_func, pars_smaxda, jobname = 'TMBrun_smaxda',
                    nodes = 100, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))


res_smaxda <- get_slurm_out(sjobtmb_smaxda, outtype = 'table', wait = TRUE)



head(res_smaxda, 3)

saveRDS(res_smaxda[res_smaxda$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_smaxda1.rds")
saveRDS(res_smaxda[res_smaxda$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_smaxda2.rds")
saveRDS(res_smaxda, file = "res_smaxda.rds")

.rslurm_id <- 56
.rslurm_istart <- (.rslurm_id)* 100 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 100, nrow(pars_smaxda))
rslurm_res_smaxda56<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_smaxda56[[i-5600]]<-tmb_func(path=".",
  a=pars_smaxda$a[i],
  u=pars_smaxda$u[i])
}
result_smaxda_56<-do.call(rbind, rslurm_res_smaxda56)
saveRDS(result_smaxda_56, file = "res_smaxda_56.rds")

#============================================================================
#sensitivity sigma scenarios low
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

pars_siglow<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_siglow <- slurm_apply(tmb_func, pars_siglow, jobname = 'TMBrun_siglow',
                    nodes = 100, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_siglow <- get_slurm_out(sjobtmb_siglow, outtype = 'table', wait = TRUE)


saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_siglow1.rds")
saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_siglow2.rds")
saveRDS(res_siglow, file = "res_siglow.rds")


#The following files are missing: results_5.RDS, results_26.RDS

.rslurm_id <- 5
.rslurm_istart <- (.rslurm_id)* 90 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 90, nrow(pars_siglow))
rslurm_res_siglow5<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_siglow5[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_siglow$a[i],
  u=pars_siglow$u[i])
}
result_siglow_5<-do.call(rbind, rslurm_res_siglow5)
saveRDS(result_siglow_5, file = "res_siglow_5.rds")



.rslurm_id <- 26
.rslurm_istart <- (.rslurm_id)* 90 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 90, nrow(pars_siglow))
rslurm_res_siglow26<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_siglow26[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_siglow$a[i],
  u=pars_siglow$u[i])
}
result_siglow_26<-do.call(rbind, rslurm_res_siglow26)
saveRDS(result_siglow_26, file = "res_siglow_26.rds")


#============================================================================
#sensitivity sigma scenarios med
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

pars_sigmed<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_sigmed <- slurm_apply(tmb_func, pars_sigmed, jobname = 'TMBrun_sigmed',
                    nodes = 100, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_sigmed <- get_slurm_out(sjobtmb_sigmed, outtype = 'table', wait = TRUE)


saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_sigmed1.rds")
saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_sigmed2.rds")
saveRDS(res_sigmed, file = "res_sigmed.rds")


#The following files are missing: results_73.RDS

.rslurm_id <- 73
.rslurm_istart <- (.rslurm_id)* 90 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 90, nrow(pars_sigmed))
rslurm_res_sigmed73<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_sigmed73[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_sigmed$a[i],
  u=pars_sigmed$u[i])
}
result_sigmed_73<-do.call(rbind, rslurm_res_sigmed73)
saveRDS(result_sigmed_73, file = "res_sigmed_73.rds")



#===========================================================================================
#---------------------------------------------------------------------------------------------------------
#stan version
library(cmdstanr)
library(rslurm)
library(samEst)
source("R/stan_func.R")


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


simPars <- read.csv("data/generic/SimPars.csv")

tst <- stan_func(path=".",
  a=1,
  u=1)


#I can run a maximum of 6 scenarios at a time, otherwise R cannot run the results in. 

#base case scn 1-6
pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)/2),each=10),
  u=1:10)

sjobstan <- slurm_apply(stan_func, pars, jobname = 'stanrun',
                    nodes = 50, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst", "cmdstanr","rlang"),
                    rscript_path = "/fs/vnas_Hdfo/comda/dag004/homey/cluster-tvsimest/",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/dag004/Rlib/4.1",
                    global_objects=c("simPars", "mod1", "mod2", "mod3",
                      "mod4","mod5","mod6","mod7","mod8"))



#AFTER JOB IS DONE IMPORT  the results
resstan <- get_slurm_out(sjobstan, outtype = 'table', wait = FALSE)

saveRDS(resstan, file = "resstan.rds")
saveRDS(resstan[resstan$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resstan1.rds")
saveRDS(resstan[resstan$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "resstan2.rds")

head(resstan, 3)



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



#stan lfo runs####
library(rslurm)
source("R/stan_lfo_func.R") #stan lfo function
library(cmdstanr)

#load in cmdstanr models for LFO
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1loo.stan")
mod1lfo=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2loo.stan")
mod2lfo=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3loo.stan")
mod3lfo=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4loo.stan")
mod4lfo=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5loo.stan")
mod5lfo=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6loo.stan")
mod6lfo=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7loo.stan")
mod7lfo=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8loo.stan")
mod8lfo=cmdstanr::cmdstan_model(file8)

#simulation parameters
simPars <- read.csv("data/generic/SimPars.csv")

#test function- this takes awhile....
#tst <- stan_lfo(path=".",
#                a=1,
#                u=1)

#test some pars
pars<-data.frame(path="..",
                   a=rep(seq_len(nrow(simPars)),each=5),
                   u=1:5)
#slurm job
sjobstan_1 <- slurm_apply(stan_lfo, pars, jobname = 'stanrun1',
                            nodes = 300, cpus_per_node = 1, submit = FALSE,
                            pkgs=c("cmdstanr"),
                            rscript_path = "/gpfs/fs7/dfo/hpcmc/comda/dag004/results/cluster-tvsimest/",
                            libPaths="/gpfs/fs7/dfo/hpcmc/comda/dag004/Rlib/4.1",
                            global_objects=c("simPars", "mod1lfo", "mod2lfo", "mod3lfo",
                                             "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))



#empirical - lfo stan runs####
library(cmdstanr)
library(rslurm)
library(samEst)
source("R/emp_lfo.R")

data<- read.csv("data/emp/salmon_productivity_compilation_may2023.csv")
stocks<- read.csv("data/emp/all_stocks_info_may2023.csv")

###Load in data####
set_cmdstan_path("/fs/vnas_Hdfo/comda/dag004/.cmdstan/cmdstan-2.31.0")
file1lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m1loo.stan")
mod1lfo=cmdstanr::cmdstan_model(file1lfo)
file2lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m2loo.stan")
mod2lfo=cmdstanr::cmdstan_model(file2lfo)
file3lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m3loo.stan")
mod3lfo=cmdstanr::cmdstan_model(file3lfo)
file4lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m4loo.stan")
mod4lfo=cmdstanr::cmdstan_model(file4lfo)
file5lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m5loo.stan")
mod5lfo=cmdstanr::cmdstan_model(file5lfo)
file6lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m6loo.stan")
mod6lfo=cmdstanr::cmdstan_model(file6lfo)
file7lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m7loo.stan")
mod7lfo=cmdstanr::cmdstan_model(file7lfo)
file8lfo=file.path(cmdstanr::cmdstan_path(),'sr models', "m8loo.stan")
mod8lfo=cmdstanr::cmdstan_model(file8lfo)


#Remove stocks with less than 15 years of recruitment data
stocks_f=subset(stocks,n.years>=16) #264 stocks
stocks_f$stock.name=gsub('/','_',stocks_f$stock.name)
stocks_f$stock.name=gsub('&','and',stocks_f$stock.name)

data_f=subset(data,stock.id %in% stocks_f$stock.id)
length(unique(data_f$stock.id)) #264
stocks_f$stock.id2=seq(1:nrow(stocks_f))
data_f$stock.id2=stocks_f$stock.id2[match(data_f$stock.id,stocks_f$stock.id)]

if(any(data_f$spawners==0)){data_f$spawners=data_f$spawners+1;data_f$logR_S=log(data_f$recruits/data_f$spawners)}
if(any(data_f$recruits==0)){data_f$recruits=data_f$recruits+1;data_f$logR_S=log(data_f$recruits/data_f$spawners)}
data_f$logR_S=log(data_f$recruits/data_f$spawners)

pars=data.frame(u=1:nrow(stocks_f))

sjobstan_emp <- slurm_apply(emp_lfo, pars, jobname = 'emprun',
                         nodes = 300, cpus_per_node = 1, submit = FALSE,
                         pkgs=c("samEst", "cmdstanr"),
                         rscript_path = "/home/dag004/Documents/cluster-tvsimest",
                         libPaths="/gpfs/fs7/dfo/hpcmc/comda/dag004/Rlib/4.1",
                         global_objects=c("data_f","stocks_f", "mod1lfo", "mod2lfo", "mod3lfo",
                                          "mod4lfo","mod5lfo","mod6lfo","mod7lfo","mod8lfo"))

#log_a sensitivity

#log_a sensitivity - half smax

#Smax sensitivity

#Smax sensitivity

