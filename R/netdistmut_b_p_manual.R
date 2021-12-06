
netdistmut_b_p_manual<-function(data,glength,iterations=500){
  library(ape)
  library(adegenet)

  ## read fasta and  distance matrix

  # dna distances
  dist<-dist.dna(data, model = "raw", variance = FALSE, gamma = FALSE, pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE)

  # changing the appropriate column names for later merge
  gendist<-as.matrix(dist)
  dim(gendist)

  rownames(gendist)<-namechanger103col_b_p(rownames(gendist) )
  colnames(gendist)<-namechanger103col_b_p(colnames(gendist) )

  gendist<-gendist[ order(rownames(gendist)),order(colnames(gendist))]

  ## relative distances

  coldist<-gendist[which(rownames(gendist)=="Col"),]
  length(coldist)
  # coldist<-coldist[-11]  # remove col and its dinstance to itself
  coldist<-coldist[-which(colnames(gendist)=="Col") ]
  coldist.mat<-matrix(ncol=1,coldist)


  thejk<-grep("JK",names(coldist))
  themodern<-seq_len(length(names(coldist)))[-thejk]
  length(themodern)
  length(thejk)

  moderndist<-coldist.mat[themodern]
  hist(moderndist)
  historicdist<-coldist.mat[thejk]
  hist(historicdist)


  tempmatM<-matrix(rep(moderndist,27),nrow=27,ncol=76,byrow=T)
  tempmatH<-matrix(rep(historicdist,76),nrow=27,ncol=76,byrow=F)
  reldist<-tempmatM-tempmatH
  hist(reldist)
  dim(reldist)


  ## Get the time distances

  reltime<-net_timedist_b_p()

  plot(reldist~reltime)

  ## Get the regression and mutation

  dataregression<-data.frame(reldist=as.numeric(reldist),reltime=as.numeric(reltime))

  ## send the bootstrap
  lengthpoints<-dim(dataregression)[1]
  slopes<-c()

  for (i in 1:iterations){
    sampledpoints<-sample(x = lengthpoints,size = lengthpoints,replace=T)
    newdataregression=dataregression[sampledpoints,]
    lmmodel<-lm(newdataregression$reldist~newdataregression$reltime)
    slope<-lmmodel$coefficients[2]
    slopes<-append(x = slopes,values = slope)
  }
  slopes<-slopes*(dim(data)[2]/glength)

  mslope<-mean(slopes)
  sdslope<-sd(slopes)
  confslope<-quantile(slopes,probs=c(0.025,0.975))

  return(matrix(c(mslope,sdslope,confslope[1],confslope[2] ) ))

}
