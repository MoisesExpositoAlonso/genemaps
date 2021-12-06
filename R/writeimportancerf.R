writeimportance.randomForest <- function(myrf, filename,path='tables'){
require(randomForest)  

    sortedimportance<-data.frame(importance(myrf) )
    sortedimportance<-rbind(sortedimportance,median(myrf$rsq))
    # sortedimportance$r2<-tail(myrf$rsq)
    sortedimportance$variable<-c(row.names(sortedimportance)[-nrow(sortedimportance)], 'R2')
    moiR::write.tsv(sortedimportance,file = paste0(path, '/',filename,"randomforestimportance.tsv"))

  }