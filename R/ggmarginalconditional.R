ggmarginalconditional<-function(x,y,ylab='Conditional effects', xlab='Marginal effects',colour='darkgrey',se=F)
{

d<-data.frame(x,y)
d$xnozero<- x ;
d$xnozero[d$xnozero==0]<-NA
d$ynozero<- y;
d$ynozero[d$ynozero==0]<-NA

head(d)

tmp<-
  ggplot(data=d)+
  geom_point(aes(x=x, y=y),shape=19) +
  xlab(xlab)+ylab(ylab) +
  stat_smooth(data=d,method="glm",se=se,colour=colour,aes(x=xnozero,y=ynozero)) +
  annotate("text",  x=Inf, y = Inf, label = lm_eq(d$ynozero,d$xnozero), vjust=1, hjust=1,parse=TRUE)

return(tmp)

}

