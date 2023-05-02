library(rslurm)
library(samEst)
source("R/tmb_func.R")



my_get_slurm_out <- function (slr_job_name, nodes.list=0:50, outtype = "raw") 
{
       
    res_files <- paste0("results_", nodes.list, ".RDS")
    tmpdir <- paste0("_rslurm_", slr_job_name)
    missing_files <- setdiff(res_files, dir(path = tmpdir))
    if (length(missing_files) > 0) {
        missing_list <- paste(missing_files, collapse = ", ")
        warning(paste("The following files are missing:", missing_list))
    }
    res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
    if (length(res_files) == 0) 
        return(NA)
    
        slurm_out <- lapply(res_files, readRDS)
    
    slurm_out <- do.call(c, slurm_out)
    if (outtype == "table") {
        slurm_out <- as.data.frame(do.call(rbind, slurm_out))
    }
    return(slurm_out)
}



#============================================================================
#sbase case scenarios

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



#run 1

saveRDS(res_a, file = "res_sensitivity_a.rds")



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


simPars <- read.csv("data/sensitivity_halfSmax/SimPars.csv")

pars_asmax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_asmax <- slurm_apply(tmb_func, pars_asmax, jobname = 'TMBrun_asmax',
                    nodes = 50, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))





#---------------------------------------------------------------------------------------------------------
#stan version


library(cmdstanr)
library(rslurm)
library(samEst)


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


simPars <- read.csv("data/generic/SimPars.csv")




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


#log_a sensitivity

#log_a sensitivity - half smax

#Smax sensitivity

#Smax sensitivity

