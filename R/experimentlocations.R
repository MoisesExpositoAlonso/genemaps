experimentlocations<-function(which=NULL){
m<-c('lat'=40.408057,'lon'= -3.835382)
t<-c('lat'=48.545860,'lon'= 9.042488)

if(is.null(which)) return(cbind(m,t))
else if(which=='m') return(m)
else if(which =='t') return(t)
else stop('Do not know the abbreviature')
}