library(ggplot2)
library(gridExtra)
library(dplyr)
source("R/utils.R")

mytheme = list(
  theme_classic(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="right", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

simPar <- read.csv("data/generic/SimPars.csv")

res<-readRDS(file = "outs/res_retro.rds")
res$parameter[res$parameter=="Smax"]<-"smax"
res<- res[res$conv_warning==0,]

resparam<-res[res$parameter%in%c("logalpha","smax","sigma","smsy","sgen","umsy"),]
resparam<-resparam[resparam$model%in%c("hmma","hmmb","rwa","rwb"),]
resparam$parameter=factor(resparam$parameter)
resparam$parameter=droplevels(resparam$parameter)
resparam$model=factor(resparam$model)
resparam$model=droplevels(resparam$model)

#resparam<- resparam[resparam$conv_warning==0,]
res2<- resparam %>% group_by(parameter,by,endyr,model) %>% summarize(median.est=median(mode,na.rm=T),q10.est=quantile(mode,0.1,na.rm=T),q90.est=quantile(mode,0.9,na.rm=T))



logalpha_retro<-res2[res2$parameter=="logalpha",]
logalpha_retro<- logalpha_retro[order(logalpha_retro$model,logalpha_retro$endyr),]

logalpha_sim=data.frame(sim=resparam$sim[match(seq(min(resparam$by),max(resparam$by)),resparam$by)],by=seq(min(resparam$by),max(resparam$by))-50)

logalpha_rwa<- logalpha_retro[logalpha_retro$model=='rwa',]

grwa<-ggplot(logalpha_rwa) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=logalpha_sim,aes(x=by,y=sim),linewidth=1.5)+
  scale_color_viridis_d("type of trend:",begin=.1, end=.8) +
  scale_fill_viridis_d("type of trend:",begin=.1, end=.8) +
  mytheme+ 
  ylab("Median estimate") +
  xlab("year")
  
grwa  

logalpha_hmma<- logalpha_retro[logalpha_retro$model=='hmma',]

ghmma<-ggplot(logalpha_hmma) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=logalpha_sim,aes(x=by,y=sim),linewidth=1.5)+
  scale_color_viridis_d("type of trend:",begin=.1, end=.8) +
  scale_fill_viridis_d("type of trend:",begin=.1, end=.8) +
  mytheme+ 
  ylab("Median estimate") +
  xlab("year")

ghmma

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  grwa
)
comb<-cowplot::plot_grid(
  grwa+ theme(legend.position="none"), ghmma+ theme(legend.position="none"),
  align = "h", axis = "bt",
  rel_widths = c(.5,.5)
)
comb2<-cowplot::plot_grid(comb, legend,
  align = "h", axis = "bt",
  rel_widths = c(.8,.2)
)
comb2
ggsave("outputs/figs/retrospective_logalpha.png",plot=stclass,width=10,height=8)

