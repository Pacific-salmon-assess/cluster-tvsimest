

#============================================================================
#sbase case scenarios


library(rslurm)
library(samEst)
source("R/tmb_func_prior_comp.R")

simPars <- read.csv("data/generic/SimPars.csv")

#base case 
#tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
#  a=5,
#  u=1)
  
tst1<-tmb_func_prior_comp(path=".",
  a=6,
  u=696)

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func_prior_comp, pars, jobname = 'TMBrunprior',
                    nodes = 240, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.3",
                    global_objects=c("simPars"))


res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)





#AFTER JOB IS DONE IMPORT  the results

saveRDS(res[res$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resbase1.rds")
saveRDS(res[res$scenario%in%simPars$scenario[floor(nrow(simPars)/2+1):nrow(simPars)],], file = "resbase2.rds")
saveRDS(res, file = "resbase.rds")



#============================================================================


# run things locally
outmats<-list()


for(i in nrow(pars)){
  a<-pars$a[i]
  u<-pars$u[i]

  outmats[[i]]<-tmb_func_prior_comp(path=".", a=a, u=u)
}


outmatsdf<- do.call(rbind,outmats)

saveRDS(outmatsdf, file = "outs/local_prior_comp.rds")