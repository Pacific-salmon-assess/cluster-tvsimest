

#============================================================================
#sbase case scenarios

library(rslurm)
library(samEst)
source("R/tmb_func_rwcompare.R")


simPars <- read.csv("data/generic/SimPars.csv")

#base case 
#tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
#  a=5,
#  u=1)
  
tst2<-tmb_func_rw_comp(path=".",
  a=3,
  u=141)
tst2

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func_rw_comp, pars, jobname = 'TMBrwcomp',
                    nodes = 250, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)
#res <- my_get_slurm_out(slr_job_name='TMBrun', nodes.list=0:39, outtype = 'table')



#AFTER JOB IS DONE IMPORT  the results

saveRDS(res[res$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "rwcompare1.rds")
saveRDS(res[res$scenario%in%simPars$scenario[floor(nrow(simPars)/2+1):nrow(simPars)],], file = "rwcompare1.rds")
saveRDS(res, file = "rwcompare.rds")

