
mode<-function(x,minS=min(x),maxS=max(x)){

  d=density(x, n=1000000, from=minS,to=maxS )

  return(d$x[which.max(d$y)])

}