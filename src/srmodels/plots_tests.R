#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================




library(ggplot2)
library(gridExtra)
source("code/utils.R")
source("code/cluster_func_plots.R")

mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res<-readRDS(file = "outs/simest/generic/res.rds")

#
#resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
#resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
#resstan<-rbind(resstan1,resstan2)
resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(res,resstan)

res$parameter[res$parameter=="smax"]<-"Smax"

resparam<-res[res$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]
resparam$bias<-(resparam$pbias/100)*resparam$sim



df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("est","sim"))

unique(df$scenario)


df_alpha_sim<- df[df$parameter=="alpha"&df$variable=="sim",]


df_alpha_est<- df[df$parameter=="alpha"&df$variable=="est",]

summarydf_alpha<-aggregate(df_alpha_est$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha<-do.call(data.frame, summarydf_alpha)

summarydf_alpha_sim<-aggregate(df_alpha_sim$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ),
    function(x) {unique(x)})


summarydf_alpha<-aggregate(df_alpha_est$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha<-do.call(data.frame, summarydf_alpha)


summarydf_alpha$scenario<-factor(summarydf_alpha$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))


summarydf_alpha_sim$scenario<-factor(summarydf_alpha_sim$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))


summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary" ),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c(       "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,3))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_alpha.png")
#MLE estimates are less biased and higher than MCMC


summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,3))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_alpha2.png")
#MLE estimates are less biased and higher than MCMC




#=======================================================
#b estimates


df_smax_sim<- df[df$parameter=="Smax"&df$variable=="sim",]


df_smax_est<- df[df$parameter=="Smax"&df$variable=="est",]

summarydf_smax<-aggregate(df_smax_est$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smax<-do.call(data.frame, summarydf_smax)

summarydf_smax_sim<-aggregate(df_smax_sim$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by ),
    function(x) {unique(x)})
head(summarydf_smax)

summarydf_smax<-aggregate(df_smax_est$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smax<-do.call(data.frame, summarydf_smax)



summarydf_smax$scenario<-factor(summarydf_smax$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))


summarydf_smax_sim$scenario<-factor(summarydf_smax_sim$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))




head(summarydf_smax)

head(summarydf_alpha_sim)
unique(summarydf_smax_sim$scenario)


summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smax)

ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(80000,300000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smax1.png")



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c("autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smax)

ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(80000,300000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smax2.png")





#=======================================================
#b estimates


df_smsy_sim<- df[df$parameter=="smsy"&df$variable=="sim",]


df_smsy_est<- df[df$parameter=="smsy"&df$variable=="est",]

summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {unique(x)})


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)



summarydf_smsy$scenario<-factor(summarydf_smsy$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))


summarydf_smsy_sim$scenario<-factor(summarydf_smsy_sim$scenario,levels=c("stationary",
                                                                   "autocorr",
                                                                   "sigmaShift",
                                                                   "decLinearProd",
                                                                   "sineProd", 
                                                                   "regimeProd", 
                                                                   "shiftProd",                       
                                                                   "decLinearCap",
                                                                   "regimeCap",
                                                                   "shiftCap",                      
                                                                   "regimeProdCap",         
                                                                   "decLinearProdshiftCap"  ))






summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smsy)

ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(80000,300000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smsy1.png")



summarydf_smsy_sim2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smsy2<-summarydf_smsy[summarydf_smsy$scenario%in%c("autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smsy)

ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method))+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(80000,300000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smsy2.png")





