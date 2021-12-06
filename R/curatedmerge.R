curatedmerge<-function(genomes,Y,name='Y'){
  ids<-genomes$fam$sample.ID
  ids<-fc(ids)
  Y[,1]<-fc(Y[,1])


  Yvar<-t(as.matrix(sapply(ids,function(i) Y[fn(which(Y[,1] == i)),-c(1)])))
  colnames(Yvar)<-name

  genomes$fam <- cbind(genomes$fam,Yvar)

  return(genomes)
}
