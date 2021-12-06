
testcolor<-function(mypalette=c('grey100','grey80')){
  x=seq(1,length(mypalette))
  y=rep(1,length(mypalette))
  plot(y=y,x=x,col=mypalette,pch=19, cex=5)

}

