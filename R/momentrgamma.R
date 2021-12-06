momentrgamma<-function(n,u,v){
  rgamma(n,shape = u^2/v ,rate=u/v)
}
