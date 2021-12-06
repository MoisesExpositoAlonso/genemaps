
timestart<-function(){
ptm <- proc.time()
return(ptm[1])
}

timeend<-function(thestart=timestart){
  message(paste('time in sec', c(proc.time()[1] - thestart),
            '\ntime in min', c(proc.time()[1] - thestart)/60,'\n'  ))
}
