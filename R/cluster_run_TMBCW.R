

#============================================================================
#sbase case scenarios

library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/generic/SimPars.csv")

#base case 
#tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
#  a=5,
#  u=1)
  
 

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)




sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)
#res <- my_get_slurm_out(slr_job_name='TMBrun', nodes.list=0:39, outtype = 'table')



,

#AFTER JOB IS DONE IMPORT  the results

saveRDS(res$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)], file = "resbase1.rds")
saveRDS(res$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)], file = "resbase2.rds")
saveRDS(res, file = "resbase.rds")


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
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res_a <- get_slurm_out(sjobtmb_a, outtype = 'table', wait = TRUE)

head(res_a, 3)


seq1<-seq_len(round(nrow(simPars)/4))
seq2<-(round(nrow(simPars)/4)+1):round(nrow(simPars)/2)
seq3<-(round(nrow(simPars)/2)+1):(round(nrow(simPars)/2)+round(nrow(simPars)/4))
seq4<-(round(nrow(simPars)/2)+round(nrow(simPars)/4)+1):nrow(simPars)



saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_a1.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_a2.rds")

saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq1],], file = "res_aq1.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq2],], file = "res_aq2.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq3],], file = "res_aq3.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq4],], file = "res_aq4.rds")



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
#sensitivity a  with double alpha scenarios
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sensitivity/SimPars_da.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_a <- slurm_apply(tmb_func, pars_a, jobname = 'TMBrun_a_da',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res_ada <- get_slurm_out(sjobtmb_a, outtype = 'table', wait = TRUE)

head(res_ada, 3)

saveRDS(res_ada[res_a$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_ada1.rds")
saveRDS(res_ada[res_a$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_ada2.rds")
saveRDS(res_ada, file = "res_ada.rds")



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
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/tvsimest/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))


res_smax <- get_slurm_out(sjobtmb_smax, outtype = 'table', wait = TRUE)


saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_smax1.rds")
saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[(nrow(simPars)/2+1):nrow(simPars)],], file = "res_smax2.rds")
saveRDS(res_smax, file = "res_smax.rds")

#run 95 that failed
.rslurm_id <- 240
.rslurm_istart <- (.rslurm_id)* 40 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 40, nrow(pars_smax))
rslurm_res_siglow240<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_siglow240[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_smax$a[i],
  u=pars_smax$u[i])
}
result_siglow_240<-do.call(rbind, rslurm_res_siglow105)
saveRDS(result_smax_240, file = "res_smax_240.rds")
saveRDS(result_smax_240, file = "C:/Users/worc/Documents/timevar/simest-tv/outs/simest/sSmax_sensitivity/res_smax_240.rds")

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
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_siglow <- get_slurm_out(sjobtmb_siglow, outtype = 'table', wait = TRUE)

round(nrow(simPars)/2)
round(nrow(simPars)/2)+1

saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_siglow1.rds")
saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_siglow2.rds")
saveRDS(res_siglow, file = "res_siglow.rds")


#The following files are missing: results_5.RDS, results_26.RDS

#run 95 that failed
.rslurm_id <- 105
.rslurm_istart <- (.rslurm_id)* 36 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 36, nrow(pars_siglow))
rslurm_res_siglow105<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_siglow105[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_siglow$a[i],
  u=pars_siglow$u[i])
}
result_siglow_105<-do.call(rbind, rslurm_res_siglow105)
saveRDS(result_siglow_105, file = "res_siglow_105.rds")
saveRDS(result_siglow_105, file = "C:/Users/worc/Documents/timevar/simest-tv/outs/simest/sigmalow_sensitivity/res_siglow_105.rds")





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
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_sigmed <- get_slurm_out(sjobtmb_sigmed, outtype = 'table', wait = TRUE)


saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_sigmed1.rds")
saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_sigmed2.rds")
saveRDS(res_sigmed, file = "res_sigmed.rds")


#The following files are missing: results_73.RDS
.rslurm_id <- 227
.rslurm_istart <- (.rslurm_id)* 36 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 36, nrow(pars_sigmed))
rslurm_res_sigmed227<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_sigmed227[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_sigmed$a[i],
  u=pars_sigmed$u[i])
}
result_sigmed_227<-do.call(rbind, rslurm_res_sigmed227)
saveRDS(result_sigmed_227, file = "res_sigmed_227.rds")
saveRDS(result_sigmed_227, file = "C:/Users/worc/Documents/timevar/simest-tv/outs/simest/sigmamed_sensitivity/res_sigmed_227.rds")



#============================================================================
#ER scenarios
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/genericER/SimPars_ER.csv")

pars_ER<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_er <- slurm_apply(tmb_func, pars_ER, jobname = 'TMBrun_ER',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_er <- get_slurm_out(sjobtmb_er, outtype = 'table', wait = TRUE)

seq1<-seq_len(round(nrow(simPars)/4))
seq2<-(round(nrow(simPars)/4)+1):round(nrow(simPars)/2)
seq3<-(round(nrow(simPars)/2)+1):(round(nrow(simPars)/2)+round(nrow(simPars)/4))
seq4<-(round(nrow(simPars)/2)+round(nrow(simPars)/4)+1):nrow(simPars)


saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_er1.rds")
saveRDS(res_er[res_er$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_er2.rds")

saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq1],], file = "res_erq1.rds")
saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq2],], file = "res_erq2.rds")
saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq3],], file = "res_erq3.rds")
saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq4],], file = "res_erq4.rds")

saveRDS(res_er, file = "res_er.rds")



.rslurm_id <- 138
.rslurm_istart <- (.rslurm_id)* 48 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 48, nrow(pars_ER))
rslurm_res_ER138<-list()
for(i in (.rslurm_istart):(.rslurm_iend)){
   rslurm_res_ER138[[i-(.rslurm_istart-1)]]<-tmb_func(path=".",
  a=pars_ER$a[i],
  u=pars_ER$u[i])
}
result_ER_138<-do.call(rbind, rslurm_res_ER138)
saveRDS(result_ER_138, file = "C:/Users/worc/Documents/timevar/simest-tv/outs/simest/genericER/result_ER_138.rds")
