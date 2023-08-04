check_stan_conv<- function(stansum){
  #stansum=resf4

  rhat <- sum((stansum$rhat-1)>.01)
  esstail <- stansum$ess_tail<400
  essbulk <- stansum$ess_bulk<400



  conv_mat <- data.frame(variable=stansum$variable,
    rhat=as.numeric((stansum$rhat-1)>.01),
    esstail=as.numeric(esstail),
    essbulk=as.numeric(essbulk)    ) 


    conv_mat$sumconv=apply(conv_mat[,-1],1,function(x){sum(x,na.rm=T)})

   


  return(list(convergence=sum(rhat,sum(esstail),sum(essbulk),na.rm=T),
    conv_mat=conv_mat
    ))



}


check_tmbstan_conv<- function(stansum){
  #stansum=resf4

  rhat <- sum((stansum[,colnames(stansum)=="Rhat"]-1)>.01)
  ess <- stansum[,colnames(stansum)=="n_eff"]<400
  

  conv_mat <- data.frame(variable=rownames(stansum),
    rhat=as.numeric((stansum[,colnames(stansum)=="Rhat"]-1)>.01),
    ess=as.numeric(ess),
    sumconv=as.numeric((stansum[,colnames(stansum)=="Rhat"]-1)>.01)+as.numeric(ess)
    ) 

  return(list(convergence=sum(rhat,sum(ess)),
    conv_mat=conv_mat
    ))



}
  