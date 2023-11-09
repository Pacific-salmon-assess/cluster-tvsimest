log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

#' mean of log-sum-exp
#'
#' @param x a vector 
#' @export
#' 
#' @returns mean of exponentials
#' 
#' 
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}



stan_refit<- function(mod,newdata,oos,K=2,dirichlet_prior=NULL){
  #mod = model file name - eg. 'ricker_linear_oos.stan'
  #newdata = data to train model
  #oosdata = data to predict onto
  #regime = TRUE or FALSE for regime shift models (have different data inputs)
  #K = number of potential regimes (2 or 3)
  
  if(is.null(dirichlet_prior)){
    dirichlet_prior<-matrix(1,nrow=K,ncol=K)
  }else if(nrow(dirichlet_prior)!=K |ncol(dirichlet_prior)!=K){
    stop("dirichlet_prior should be a K x K matrix")
  }
  
  oosdata=newdata[oos,]
  newdata=newdata[-oos,]
  
  dl=list(
    by=newdata$by,
    N=nrow(newdata),
    L=max(newdata$by)-min(newdata$by)+1,
    ii=newdata$by-min(newdata$by)+1,
    R_S =newdata$logRS,
    S=newdata$S,
    y_oos=oosdata$logRS,
    x_oos=oosdata$S,
    K=2,
    alpha_dirichlet= dirichlet_prior
  )
  
  r = mod$sample(data=dl,
                 seed=123,
                 chains=6,
                 iter_warmup=200,
                 iter_sampling=500,
                 refresh=0,
                 adapt_delta=0.95,
                 max_treedepth=15)
  return(r)
}


#' stan_lfo_cv function

stan_lfo_cv=function(mod,type=c('static','tv','regime'),df,L=10,K=2,dirichlet_prior=NULL){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (default 10)
  # K = number of regimes
  
  loglik_exact <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  if(type=='regime'){
    loglik_exact_1bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
    loglik_exact_3bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
    loglik_exact_5bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  }
  for (i in L:(nrow(df) - 1)){
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    if(type=='static'){
      fit_past<- stan_refit(mod=mod, newdata=df_oos, oos=i+1)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos'),format='draws_matrix'))
      loglik_exact[,i+1]<- ll$log_lik_oos
    }
    if(type=='tv'){
      fit_past<- stan_refit(mod=mod,newdata=df_oos,oos=i+1)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'),format='draws_matrix'))
      
      loglik_exact_1b[, i + 1] <-ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <-ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <-ll$log_lik_oos_5b
    }
    if(type=='regime'){
      fit_past<- stan_refit(mod=mod,newdata=df_oos,oos=i+1, K=K, dirichlet_prior=dirichlet_prior)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'),format='draws_matrix'))
      loglik_exact_1b[, i + 1] <- ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <- ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <- ll$log_lik_oos_5b
    }
  }
  
  if(type=='static'){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    r=exact_elpds
  }
  if(type=='tv'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  if(type=='regime'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  return(r)
}


emp_lfo<- function(u){
  #data_f
  dat <- data_f[data_f$stock.id2==unique(data_f$stock.id2)[u],]
  dat <- dat[complete.cases(dat$spawners),]
  df <- data.frame(by=dat$broodyear,
                   S=dat$spawners,
                   R=dat$recruits,
                   logRS=log(dat$recruits/dat$spawners),
                   pSmax_mean=max(dfdat$spawners)*0.5,
                   pSmax_sig=max(dfdat$spawners)*0.5)
  
  
  options(mc.cores = 6)
  
  #LFO cross-validation
  #model 1 - static Ricker
  lfostatic<-stan_lfo_cv(mod=mod1lfo,type='static',df=df,L=10)
  #model 2 - static autocorrelated Ricker
  lfoac<- stan_lfo_cv(mod=mod2lfo,type='static',df=df,L=10)
  #model 3 - dynamic productivity Ricker
  lfoalpha<-stan_lfo_cv(mod=mod3lfo,type='tv',df=df,L=10)
  #model 4 - dynamic capacity Ricker
  lfobeta<- stan_lfo_cv(mod=mod4lfo,type='tv',df=df,L=10)
  #model 5 - dynamic productivity & capacity Ricker
  lfoalphabeta<- stan_lfo_cv(mod=mod5lfo,type='tv',df=df,L=10)
  #model 6 - productivity regime shift - 2 regimes
  lfohmma<- stan_lfo_cv(mod=mod6lfo,type='regime',df=df,L=10,K=2,dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2))
  #model 7 - capacity regime shift
  lfohmmb<- stan_lfo_cv(mod=mod7lfo,type='regime',df=df,L=10,K=2,dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2))
  #model 8 - productivity and capacity regime shift
  lfohmm<- stan_lfo_cv(mod=mod8lfo,type='regime',df=df,L=10,K=2,dirichlet_prior=matrix(c(2,1,1,2),ncol=2,nrow=2))
  
  dflfo<- data.frame(parameter="LFO",
                     stock.id=u,
                     stock=unique(data_f$stock.name)[u] ,
                     model=c("simple",
                             "autocorr",
                             "rwa","rwb","rwab",
                             "rwa_last3","rwb_last3","rwab_last3",
                             "rwa_last5","rwb_last5","rwab_last5",
                             "hmma", "hmmb","hmmab",
                             "hmma_last3", "hmmb_last3","hmmab_last3",
                             "hmma_last5", "hmmb_last5","hmmab_last5"),
                     est=c(sum(lfostatic), 
                           sum(lfoac), 
                           sum(lfoalpha[1,]), 
                           sum(lfoalpha[2,]), 
                           sum(lfoalpha[3,]), 
                           sum(lfobeta[1,]), 
                           sum(lfobeta[2,]), 
                           sum(lfobeta[3,]), 
                           sum(lfoalphabeta[1,]), 
                           sum(lfoalphabeta[2,]), 
                           sum(lfoalphabeta[3,]),    
                           sum(lfohmma[1,]),
                           sum(lfohmma[2,]),
                           sum(lfohmma[3,]),
                           sum(lfohmmb[1,]),
                           sum(lfohmmb[2,]),
                           sum(lfohmmb[3,]),
                           sum(lfohmm[1,]),
                           sum(lfohmm[2,]),
                           sum(lfohmm[3,])
                     ))

  return(dflfo)
  
}