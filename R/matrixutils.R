

seemat<-function(m,howmuch=10){
  m[1:howmuch, 1:howmuch]
}

subsetmat<-function(m,cols){

  tosub<-which(colnames(m) %in% cols)

  stopifnot(length(tosub)>0)

  m[tosub,tosub]
}


mh<-function(x,start=1,end=5){
  x[start:end,start:end]

}

