
#============================================================================
#all scenarios

library(rslurm)
library(samEst)
source("R/tmb_func.R")

simPars1 <- read.csv("data/generic/SimPars.csv")
simPars2 <- read.csv("data/sensitivity/SimPars.csv")
simPars3 <- read.csv("data/sensitivity_halfSmax/SimPars.csv")
simPars4 <- read.csv("data/Smax_sensitivity/SimPars.csv")
simPars5 <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")
simPars6 <- read.csv("data/sigmalow_sensitivity/SimPars.csv")
simPars7 <- read.csv("data/sigmamed_sensitivity/SimPars.csv")
simPars8 <- read.csv("data/genericER/SimPars_ER.csv")
simPars9 <- read.csv("data/generic_biascorr/SimPars_biascorr.csv")


simParsall<-rbind(simPars1,simPars2,simPars3,simPars4,simPars5,simPars6,simPars7,simPars8,simPars9)



parsall<-data.frame(path="..",
  a=rep(seq_len(nrow(simParsall)),each=1000),
  u=1:1000)


sjobtmball <- slurm_apply(tmb_func, parsall, jobname = 'TMBrunall',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simParsall"))

resall <- get_slurm_out(sjobtmball, outtype = 'table', wait = TRUE)

#rowsby scenario
nrow(simParsall)

#base
saveRDS(resall[resall$scenario%in%simPars1$scenario,], file = "resbase.rds")
saveRDS(resall[resall$scenario%in%simPars1$scenario[seq_len(nrow(simPars1)/2)],],], file = "resbase1.rds")
saveRDS(resall[resall$scenario%in%simPars1$scenario[(nrow(simPars1)/2+1):nrow(simPars1)],],], file = "resbase2.rds")

#sensitivity a
saveRDS(resall[resall$scenario%in%simPars2$scenario,], file = "res_a.rds")
saveRDS(resall[resall$scenario%in%simPars2$scenario[seq_len(nrow(simPars2)/2)],], file = "res_a1.rds")
saveRDS(resall[resall$scenario%in%simPars2$scenario[(nrow(simPars2)/2+1):nrow(simPars2)],], file = "res_a2.rds")

#sensitivity a scenarios half smax
saveRDS(resall[resall$scenario%in%simPars3$scenario,], file = "res_asmax.rds")
saveRDS(resall[resall$scenario%in%simPars3$scenario[seq_len(nrow(simPars3)/2)],], file = "res_asmax1.rds")
saveRDS(resall[resall$scenario%in%simPars3$scenario[(nrow(simPars3)/2+1):nrow(simPars3)],], file = "res_asmax2.rds")

#smax scenarios 
saveRDS(resall[resall$scenario%in%simPars4$scenario,], file = "res_smax.rds")
saveRDS(resall[resall$scenario%in%simPars4$scenario[seq_len(nrow(simPars4)/2)],], file = "res_smax1.rds")
saveRDS(resall[resall$scenario%in%simPars4$scenario[(nrow(simPars4)/2+1):nrow(simPars4)],], file = "res_smax2.rds")

#smax scenarios double alpha
saveRDS(resall[resall$scenario%in%simPars5$scenario,], file = "res_smaxda.rds")
saveRDS(resall[resall$scenario%in%simPars5$scenario[seq_len(nrow(simPars5)/2)],], file = "res_smaxda1.rds")
saveRDS(resall[resall$scenario%in%simPars5$scenario[(nrow(simPars5)/2+1):nrow(simPars5)],], file = "res_smaxda2.rds")

#sensitivity sigma scenarios low
saveRDS(resall[resall$scenario%in%simPars6$scenario,], file = "res_siglow.rds")
saveRDS(resall[resall$scenario%in%simPars6$scenario[seq_len(nrow(simPars6)/2)],], file = "res_siglow1.rds")
saveRDS(resall[resall$scenario%in%simPars6$scenario[(nrow(simPars6)/2+1):nrow(simPars6)],], file = "res_siglow2.rds")

#sensitivity sigma scenarios med
saveRDS(resall[resall$scenario%in%simPars7$scenario,], file = "res_sigmed.rds")
saveRDS(resall[resall$scenario%in%simPars7$scenario[seq_len(nrow(simPars7)/2)],], file = "res_sigmed1.rds")
saveRDS(resall[resall$scenario%in%simPars7$scenario[(nrow(simPars7)/2+1):nrow(simPars7)],], file = "res_sigmed2.rds")


#sensitivity baseER
saveRDS(resall[resall$scenario%in%simPars8$scenario,], file = "resbase_ER.rds")
saveRDS(resall[resall$scenario%in%simPars8$scenario[seq_len(nrow(simPars8)/2)],], file = "resbase_ER1.rds")
saveRDS(resall[resall$scenario%in%simPars8$scenario[(nrow(simPars8)/2+1):nrow(simPars8)],], file = "resbase_ER2.rds")


#sensitivity base_biascorr
saveRDS(resall[resall$scenario%in%simPars9$scenario,], file = "resbase_biascorr.rds")
saveRDS(resall[resall$scenario%in%simPars9$scenario[seq_len(nrow(simPars9)/2)],], file = "resbase_biascorr1.rds")
saveRDS(resall[resall$scenario%in%simPars9$scenario[(nrow(simPars9)/2+1):nrow(simPars9)],], file = "resbase_biascorr2.rds")

#============================================================================
#sbase case scenarios

library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/generic/SimPars.csv")


#base case 
tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
  a=5,
  u=1)
  
 

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)




sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 200, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
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
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
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
#run 58103rm_id <- 1
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
