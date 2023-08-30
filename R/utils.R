
postmode<-function(x,minS=min(x),maxS=min(max(x),median(x)*20)){

  d=density(x, n=1000000, from=minS,to=maxS )

  return(d$x[which.max(d$y)])

}