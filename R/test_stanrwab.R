#==========================================================
#Script tro ttest rwab sstan based models
#Catarina wor 
#==========================================================

library(ggplot2)
library(cmdstanr)
library(rstan)
library(samEst)
source("R/stan_func.R")


mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



file5=file.path(cmdstanr::cmdstan_path(),'srmodels', "m5f.stan")
mod5=cmdstanr::cmdstan_model(file5)


#---------------------------------------------------------------------------------------------------
#stan  base scenarios


simPars <- read.csv("data/generic/SimPars.csv")

path="."
a=4
u=19

  
  
allsimest <- list()
  simData <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2)
  )


  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))
  
 
  
rstan_options(auto_write = TRUE)
options(mc.cores = (parallel::detectCores()-1))


f5 <- mod5$sample(data=df,
                    seed = 123,
                    chains = 30, 
                    iter_warmup = 4000,
                    iter_sampling = 30000,
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)
  




    #====================================================
    #compare rw a between, all forms of model fit


file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ja.stan")
mod3=cmdstanr::cmdstan_model(file3)



file3noadj=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f.stan")
mod3noadj=cmdstanr::cmdstan_model(file3noadj)




simPars <- read.csv("data/generic/SimPars.csv")

path="."
a=4
ulist=19:39

dfal<-list()

for(i in seq_along(ulist)){

u=ulist[i]

 
allsimest <- list()
  simData <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2)
  )

  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))



rwa_tmb <-ricker_rw_TMB(data=df_tmb ,tv.par="a",sig_p_sd=1)

kfa_tmb <- ricker_kf_TMB(data=df_tmb )

f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6,#30 
                    iter_warmup = 500, # 4000
                    iter_sampling = 5000, #30000
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)

 stanlog_a<-f3$summary(variables=c('log_a'),'median')$median




tmbstan<-ricker_rw_TMBstan(data=df_tmb, chains=6, #30,
  iter=2000, #30000
  tv.par="a",
   warmup = 500) # 4000



dfal[[i]]<-data.frame(param="a",
  model=rep(c("rwa_tmb","kf_tmb","stan","rwa_tmbstan","sim"),each=40),
  by=1:40,
  iteration=u,
  scenario=simPars$scenario[a],
  value=c(rwa_tmb$alpha,
    kfa_tmb$alpha,
    f3$summary(variables=c('log_a'),'mean')$mean,
    apply(tmbstan$alpha,1,mean),
    dat$alpha))
    



}  
 
dfa<-do.call(dfal,rbind)
f3$summary(variables=c('log_a'),)
  
  



 
#not currently working because I cannot extract posterior estimates of the laplace approximation
tmbstanlaplace<-ricker_rw_TMBstan(data=df_tmb,
 chains=6,
 iter=10000,
 laplace=TRUE,
 tv.par="a")




dfa<-data.frame(param="a",
  model=rep(c("rwa_tmb","kf_tmb","stan","rwa_tmbstan","sim"),each=40),
  by=1:40,
  value=c(rwa_tmb$alpha,
    kfa_tmb$alpha,
    f3$summary(variables=c('log_a'),'mean')$mean,
    apply(tmbstan$alpha,1,mean),
    dat$alpha))
    #apply(tmbstanlaplace$alpha,1,median)))



ggplot(dfa)+
geom_line(aes(x=by,y=value,color=model),linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
mytheme




library(tmbstan)
runExample("simple")

init.fn <- function()
  list(beta=rnorm(2), logsdu=runif(1,0,10), logsd0=runif(1,0,1))
fit <- tmbstan(obj, chains=cores, open_progress=FALSE,
               init=init.fn, laplace=TRUE)

## There are no posterior samples for the random effects because they are
## integrated out by the LA. See Monnahan and Kristensen (2019) for discussion
## of why this would be worth doing. Typically it will be slower and less
## accurate than laplace=FALSE (the default).
names(as.data.frame(fit))





file3=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_ja.stan")
mod3=cmdstanr::cmdstan_model(file3)



file3_Smax=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f_Smax.stan")
mod3_Smax=cmdstanr::cmdstan_model(file3_Smax)


file3noadj=file.path(cmdstanr::cmdstan_path(),'srmodels', "m3f.stan")
mod3noadj=cmdstanr::cmdstan_model(file3noadj)



simPars <- read.csv("data/generic/SimPars.csv")

path="."
a=4
ulist=19:39

dfal<-list()
i=1
for(i in seq_along(ulist)){

u=ulist[i]

 
allsimest <- list()
  simData <- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=.5*max(dat$obsSpawners),
             pSmax_sig=max(dat$obsSpawners)

  )

  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))



rwa_tmb <- ricker_rw_TMB(data=df_tmb ,tv.par="a",sig_p_sd=1)

kfa_tmb <- ricker_kf_TMB(data=df_tmb )

f3 <- mod3$sample(data=df,
                    seed = 123,
                    chains = 6,#30 
                    iter_warmup = 500, # 4000
                    iter_sampling = 5000, #30000
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)

f3_Smax <- mod3_Smax$sample(data=df,
                    seed = 123,
                    chains = 6,#30 
                    iter_warmup = 500, # 4000
                    iter_sampling = 5000, #30000
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)


  df$pSmax_mean=.5*max(dat$obsSpawners)
  df$pSmax_sig=.5*max(dat$obsSpawners)


f3_Smaxp2 <- mod3_Smax$sample(data=df,
                    seed = 123,
                    chains = 6,#30 
                    iter_warmup = 500, # 4000
                    iter_sampling = 5000, #30000
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)

f3noadj <- mod3noadj$sample(data=df,
                    seed = 123,
                    chains = 6,#30 
                    iter_warmup = 500, # 4000
                    iter_sampling = 5000, #30000
                    refresh = 0,
                    adapt_delta = 0.98,
                    max_treedepth = 15)

 #stanlog_a<-f3$summary(variables=c('log_a'),'median')$median



#tmbstan<-ricker_rw_TMBstan(data=df_tmb, chains=6, #30,
#  iter=2000, #30000
#  tv.par="a",
#   warmup = 500) # 4000


#
dfal[[i]]<-data.frame(param="a",
  model=rep(c("rwa_tmb","kf_tmb","stan","stan_noadj","stan_Smax","sim"),each=40),
  by=1:40,
  iteration=u,
  scenario=simPars$scenario[a],
  value=c(rwa_tmb$alpha,
    kfa_tmb$alpha,
    f3$summary(variables=c('log_a'),'median')$median,
    f3noadj$summary(variables=c('log_a'),'median')$median,
    f3_Smax$summary(variables=c('log_a'),'median')$median,
    f3_Smaxp2$summary(variables=c('log_a'),'median')$median,
    dat$alpha))
    


#ggplot(dfa)+
#geom_line(aes(x=by,y=value,color=model),linewidth=1.2)+
#scale_color_viridis_d(begin=.1, end=.8) +
#mytheme




}  
 
dfa<-do.call(rbind, dfal)
dfa<-dfa[dfa$iteration<35,]
head(dfa)
summary(dfa$iteration)

ggplot(dfa[dfa$iteration==21,])+
geom_line(aes(x=by,y=value,color=model, group=interaction(model, iteration)),linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
mytheme


summarydf_alpha<-aggregate(dfa$value,by=list(scenario=dfa$scenario, 
    model=dfa$model,
    by=dfa$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha<-do.call(data.frame, summarydf_alpha)


head(summarydf_alpha)
ggplot(summarydf_alpha)+
geom_pointrange(data=summarydf_alpha,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=model))+
geom_line(data=summarydf_alpha,aes(x=by,y=x.50.,color=model),linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option="turbo")+
mytheme






plot(seq(0,495977.3*3,10000),dnorm(seq(0,495977.3*3,10000),247988.7, 495977.3/2))
abline(v=(495977.3/2))
abline(v=(495977.3/2))


plot((seq(0,495977.3*2,10000)),dlnorm(seq(0,495977.3*2,10000),log(247988.7*2), 1))
abline(v=(495977.3))


#======================================================================
#rwb examples

file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ip.stan")
mod4=cmdstanr::cmdstan_model(file4)

file4Smax=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_smax.stan")
mod4Smax=cmdstanr::cmdstan_model(file4Smax)


file4_ja=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ja.stan")
mod4_ja=cmdstanr::cmdstan_model(file4_ja)


library(TMB)
#compile tmb models

compile("src/srmodels/Ricker_tvlogb_centered.cpp")
dyn.load("src/srmodels/Ricker_tvlogb_centered.dll")


compile("src/srmodels/Ricker_tvlogb_less1.cpp")
dyn.load("src/srmodels/Ricker_tvlogb_less1.dll")





path="."
a=5
u=4

 
simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=max(dat$obsSpawners)*.5,
             pSmax_sig=max(dat$obsSpawners)*.5
  )
  dfsip<-df
  dfsip$pSmax_sig<-max(dat$obsSpawners)*.25

  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  
  

  #prior values for TMB
  Smax_mean<-(max(df_tmb$S)*.5)
  Smax_sd<-Smax_mean
  Smax_sd_sip<-((max(df_tmb$S)*.5))*2

  logbeta_pr_sig=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_pr=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
 
  logbeta_sippr_sig=sqrt(log(1+((1/ Smax_sd_sip)*(1/ Smax_sd_sip))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_sippr=log(1/(Smax_mean))-0.5*logbeta_sippr_sig^2
 


  f4_ip <- mod4_ip$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 500,
                    iter_sampling = 5000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  resf4_ip<-f4_ip$summary()
  resf4_ip_conv<-check_stan_conv(resf4_ip)

  ptvb_ip <- ricker_rw_TMB(data=df_tmb, tv.par='b',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
  
  stvb_ip<-ricker_rw_TMBstan(data=df_tmb, tv.par='b',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig,
    chains=10,
    iter=5000 ,
    warmup =500,
     control = list(adapt_delta = 0.99))

  apply(stvb_ip$Smax,1,median)
  
length(stvb_ip$alpha)

  #centered around 0 tmb and tmbstan
  
  initlm<-lm(logRS~S, data=df_tmb)

  

  tmb_data <- list(
    obs_S = df_tmb$S,
    obs_logRS = df_tmb$logRS,
    priors_flag=1,
    stan_flag=0,
    sig_p_sd=1,
    sigb_p_sd=1,
    logb_p_mean=logbeta_pr,
    logb_p_sd=logbeta_pr_sig
  )

  
tmb_params_ce <- list(logbetao = log(1/(max(df_tmb$S)*.5)),
                   alpha   = max(initlm$coefficients[[1]],.5),                 
                   logsigobs = log(.6),
                   logsigb = log(.2),
                   bdev=rep(0, length(df_tmb$S)-1))
obj_cent<-MakeADFun(tmb_data,tmb_params_ce,DLL="Ricker_tvlogb_centered")
    newtonOption(obj_cent, smartsearch=FALSE)

opt_cent<-nlminb(obj_cent$par,obj_cent$fn,obj_cent$gr)
rep_cent<-obj_cent$report()


convcent<-get_convergence_diagnostics(TMB::sdreport(obj_cent))

)

tmb_params_le1 <- list(logbetao = log(1/(max(df_tmb$S)*.5)),
                   alpha   = max(initlm$coefficients[[1]],.5),                 
                   logsigobs = log(.6),
                   logsigb = log(.2),
                   logbeta=rep(log(1/(max(df_tmb$S)*.5)), length(df_tmb$S)-1))

obj_le1<-MakeADFun(tmb_data,tmb_params_le1,DLL="Ricker_tvlogb_less1")
    newtonOption(obj_cent, smartsearch=FALSE)

opt_le1<-nlminb(obj_le1$par,obj_le1$fn,obj_le1$gr)
rep_le1<-obj_le1$report()


convle1<-get_convergence_diagnostics(TMB::sdreport(obj_le1))

grep("S_max",stvb_ip$fit_summary$summary$variable)
head(stvb_ip$fit_summary$summary)

stvb_ip$Smax
  stvb_ip_conv<- check_tmbstan_conv(stansum=stvb_ip$fit_summary$summary)

  stvb_ip_conv$conv_mat$sumconv[grep("logbeta\\[",stvb_ip_conv$conv_mat$variable)]

  stvb_ip_conv$conv_mat$variable[grep("logbeta\\[",stvb_ip_conv$conv_mat$variable)]
  names(stvb_ip)


  dfsmax<-data.frame(param="smax",
  model=rep(c("rwb_tmb","rwb_tmb_le","rwb_tmb_ce","rwb_tmbstan","rwb_stan","sim"),each=40),
  by=1:40,
  iteration=u,
  scenario=simPars$scenario[a],
  value=c(ptvb_ip$Smax,
    rep_le1$Smax,
    rep_cent$Smax,
    apply(stvb_ip$Smax,1,median),
    resf4_ip[grep("S_max",resf4_ip$variable),"median"][[1]],   
    dat$capacity),
  convergence=c(rep(ptvb_ip$model$convergence + ptvb_ip$conv_problem,nrow(dat)),
    rep(convle1$conv_problem,nrow(dat)),
    rep(convcent$conv_problem,nrow(dat)),
    stvb_ip_conv$conv_mat$sumconv[grep("logbeta\\[",stvb_ip_conv$conv_mat$variable)],
    resf4_ip_conv$conv_mat$sumconv[grep("S_max\\[",resf4_ip_conv$conv_mat$variable)],
    rep(0,nrow(dat))
    ))


dfsmax$value[dfsmax$convergence!=0]<-NA


ggplot(dfsmax)+
geom_line(aes(x=by,y=value,color=model, group=interaction(model, iteration)),linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
mytheme



#======================================================================
#rwb examples
  options(mc.cores = 5)
file4=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f.stan")
mod4=cmdstanr::cmdstan_model(file4)


file4_ip=file.path(cmdstanr::cmdstan_path(),'srmodels', "m4f_ip.stan")
mod4_ip=cmdstanr::cmdstan_model(file4_ip)


library(TMB)
#compile tmb models

path="."
a=6
u=10

 
simData<- readRDS(paste0(path,"/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData[simData$iteration==u,]
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- list(by=dat$year,
             S=dat$obsSpawners,
             R=dat$obsRecruits,
             R_S=log(dat$obsRecruits/dat$obsSpawners),
             L=max(dat$year)-min(dat$year)+1,
             ii=as.numeric(as.factor(dat$year)),
             N=nrow(dat),
             K=2,
             alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),
             pSmax_mean=max(dat$obsSpawners)*.5,
             pSmax_sig=max(dat$obsSpawners)*.5
  )
  dfsip<-df
  dfsip$pSmax_sig<-max(dat$obsSpawners)*.25

  df_tmb <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

  
  

  #prior values for TMB
  Smax_mean<-(max(df_tmb$S)*.5)
  Smax_sd<-Smax_mean
  Smax_sd_sip<-((max(df_tmb$S)*.5))*2

  logbeta_pr_sig=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_pr=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
 
  logbeta_sippr_sig=sqrt(log(1+((1/ Smax_sd_sip)*(1/ Smax_sd_sip))/((1/Smax_mean)*(1/Smax_mean))))
  logbeta_sippr=log(1/(Smax_mean))-0.5*logbeta_sippr_sig^2
 


  f4_ip <- mod4_ip$sample(data=df,
                    seed = 123,
                    chains = 6, 
                    iter_warmup = 1000,
                    iter_sampling = 10000,
                    refresh = 0,
                    adapt_delta = 0.99,
                    max_treedepth = 15)
  resf4_ip<-f4_ip$summary()
  resf4_ip_conv<-check_stan_conv(resf4_ip)

 head(f4_ip)
 draws <- f4_ip$draws(format = "matrix")
 colnames(draws)[grep("S_max\\[",colnames(draws))]
 
 "S_max[40]"

bayesplot::mcmc_areas(
  f4_ip$draws(colnames(draws)[grep("S_max\\[",colnames(draws))]), 
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "median"
)+
theme_bw(12)


  ptvb_ip <- ricker_rw_TMB(data=df_tmb, tv.par='b',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig)
  
  stvb_ip<-ricker_rw_TMBstan(data=df_tmb, tv.par='b',logb_p_mean=logbeta_pr,logb_p_sd=logbeta_pr_sig,
    chains=10,
    iter=5000 ,
    warmup =500,
     control = list(adapt_delta = 0.99))
  
  stvb_ip_conv<- check_tmbstan_conv(stansum=stvb_ip$fit_summary$summary)


  dfsmax<-data.frame(param="smax",
  model=rep(c("rwb_tmb","rwb_tmbstan","rwb_stan","sim"),each=40),
  by=1:40,
  iteration=u,
  scenario=simPars$scenario[a],
  value=c(ptvb_ip$Smax,
    apply(stvb_ip$Smax,1,median),
    resf4_ip[grep("S_max",resf4_ip$variable),"median"][[1]],   
    dat$capacity),
  convergence=c(rep(ptvb_ip$model$convergence ptvb_ip$conv_problem,nrow(dat)),
    stvb_ip_conv$conv_mat$sumconv[grep("logbeta\\[",stvb_ip_conv$conv_mat$variable)],
    resf4_ip_conv$conv_mat$sumconv[grep("S_max\\[",resf4_ip_conv$conv_mat$variable)],
    rep(0,nrow(dat))
    ))


dfsmax$value[dfsmax$convergence!=0]<-NA


ggplot(dfsmax)+
geom_line(aes(x=by,y=value,color=model, group=interaction(model, iteration)),linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
mytheme

