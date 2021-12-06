
nucdiv<-function(alldifferences,numsites=1353386,totcompar=NULL){

    if(is.null(totcompar)) {
      # totcompar=combn(length(alldifferences),m = 2)
      totcompar=table(!is.na((alldifferences)))["TRUE"]
      totcompar=as.numeric(totcompar)
    }

  div=mean( alldifferences/(numsites*totcompar), na.rm = TRUE )

  return(div)

}

whichclosest<-function(vectordist,numberneighbors=10,k=NULL){
  distthreshold<-sort(vectordist)[numberneighbors+1]
  who<-as.numeric(which(vectordist <distthreshold & vectordist >0))
 if(!is.null(k) & length(who)<k){
 who<-as.numeric(names(sort(vectordist))[2:numberneighbors+1])

 }
   return(who)
}

neighbordiversity<-function(hammingm=ham,geographicm=geodist,numberneighbors=10, k=2){
  newdiv<-c()
  for (i in 1:dim(geographicm) [1]){
    who<-whichclosest(geographicm[i,],numberneighbors = numberneighbors,k=k)

    # print(who)
        if (k==2){
        nearbydistance<-ham[i,who]
        newdiv[i] <- nucdiv(alldifferences = sample(nearbydistance,size=1), totcompar = 1,numsites = 119e6) # be careful with more than 2 you have to find the diversity of all comparisons, not only from a focal genome
        }else if (k>2){
          samplewho<-sample(who,size=k,replace = F)
        nearbydistance<-ham[  samplewho,  samplewho]
        nearbydistance_compar<-nearbydistance[lower.tri(nearbydistance)]
        totcompar=length(nearbydistance_compar)
        newdiv[i] <- nucdiv(alldifferences = nearbydistance_compar, totcompar = totcompar,numsites = 119e6)
        }
  }
  hist(newdiv,col="black",border="white",main = paste("nei=",numberneighbors,"|","k=",k),breaks=20)
  return(newdiv)
}


div_perpop<-function(acclist,bigmatrix){

submat=bigmatrix[as.character(acclist) ,as.character(acclist)]

nucdiv(submat[lower.tri(submat)==T],totcompar = length(submat[lower.tri(submat)==T]) , numsites = 119e6)

}
