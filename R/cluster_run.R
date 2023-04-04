library(rslurm)
library(samEst)


simPars <- read.csv("data/generic/SimPars.csv")

source("R/tmb_func.R")


tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=1,
  u=1)
  
pars<-data.frame(path="..",
  a=rep(1:11,each=1000),
  u=1:1000)

res<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/cluster-tvsimest",
  a=1,
  u=1)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 4, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/caw001/Documents/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res <- get_slurm_out(sjobtmb, outtype = 'table', wait = FALSE)
head(res, 3)



saveRDS(res, file = "res.rds")