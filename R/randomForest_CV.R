#' Title
#'
#' @param mydata
#' @param var
#' @param mypredictors
#' @param cross
#' @param seed
#' @param myfilename
#' @param type
#' @param sink
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples

randomForest_CV<-function(mydata,var,mypredictors,cross=5,seed=1, type="regression",sink=T, out='tmptables' ,filename="",verbose=T){

# checking arguments
if(!type %in%  c('regression','classification')){message('provide an appropriate random forest type')}

# load llibraries
library(caret)
library(randomForest)

# generate the cross validation partitions
if(verbose) message(paste('creating',cross, 'partitions of the data ...'))
group<-createFolds(mydata[,'beta'] , k=cross,list=FALSE)

# set up several parameters
allpredicted=c()
varcol=which(colanmes(mydata)==var)

if(verbose) message(paste('running Random Forest...'))
rfmodels<-lapply(1:cross,FUN =
          function(l){
            message(paste('# ',l))

            train <- mydata[group != l,]
            test <- mydata[group == l,-varcol]

            if(type=="classification"){ mydata[,var] = factor(mydata[,var])
            }else{mydata[,var] = as.numeric(mydata[,var])}

            myformula<-formula(paste(var," ~ ", paste(mypredictors,collapse = '+')))

            fit <- train(myformula, data=mydata,method='rf', prox=TRUE )
            # fit <- randomForest(myformula, data=mydata,ntree=500,na.action=na.omit,keep.inbag = TRUE,importance=TRUE) # removed proximity

            mypred <-predict(fit,test)
            allpredicted<-c(allpredicted,mypred)

          return(fit)
} )

# combine all random forest
if(verbose) message("merging the random forest")
rfall<-randomForest::combine(rfmodels[[1]],rfmodels[[2]],rfmodels[[3]],rfmodels[[4]],rfmodels[[5]])

# Generate accuracy estimation
if(type=="classification"){
    myr2= table(ypredict == mydata[,phenoname])["TRUE"] / sum(table(ypredict == mydata[,phenoname]))
    if(verbose) message(paste("predictive power, proportion of right predictions=",myr2))
}else{
    myr2<-summary(lm(allpredicted~mydata[,phenoname]))$r.squared
    if(verbose) message(paste("predictive power, R2=",myr2))
}

# sink info
if(sink==T){
    sortedimportance<-data.frame(importance(rfall) )
    sortedimportance$r2<-myr2
    sortedimportance$var<-row.names(sortedimportance)
    write.tsv(sortedimportance,file = paste0(out, '/',var,"_",filename,"randomforest.tsv"))
}

return(rfall)
}


#
# randomForest_CV<-function(mydata,phenoname,cross=5,seed=1, myfilename="",type="regression",sink=T){
#
# # This function computes a random forest from a training dataset. To have replicable data we set seed at 1. Use something else if you want
# # If cross=0, then no cross validation is done.
#
# set.seed(seed)
# require(randomForest)
# require(reshape)
# require(dismo)
#
# group <- dismo::kfold(mydata, k=cross)
# print(paste("calculating random forest of ...",phenoname))
#
# rfmodels<-lapply(1:cross,FUN =
# function(l){
#   # train,mydata,mydata
#   train <- mydata[group != l,]
#   if(type=="classification"){
#     myformula<-formula(paste("factor(", phenoname,") ~ ."))
#   } else if(type=="regression"){ # so regression
#       #myformula<-formula(paste(phenoname," ~ ."))
#       myformula<-formula(paste('as.numeric(',phenoname,") ~ .")) # problems running regression with alleles, likely because they are characters 0 1 NA. need to put as.numeric
#       } else{ print("need to provide a valid type of random forest: regression or classification")}
#   fit <- randomForest(myformula, data=mydata,ntree=500,na.action=na.omit,keep.inbag = TRUE,importance=TRUE) # removed proximity
# return(fit)
# } )
#
#
# print("join random forest cross validations...")
# # merge the cross validated random forests
# # pasmodel<-paste0("rfmodels [[",c(1:cross),"]]", collapse = ",")
# # rfall<-evalparse(paste("forcecombine_randomforest(",pasmodel,")") )
# # rfall<-evalparse(paste("combine(",pasmodel,")") )
# rfall<-combine(rfmodels[[1]],rfmodels[[2]],rfmodels[[3]],rfmodels[[4]],rfmodels[[5]])
#
# ### evaluate the prediction of rf regression predicted and observed.
# ypredict<-predict(object = rfall,newdata=mydata)
#
#
# if(type=="classification"){
# print("as binary classification ...")
# myr2= table(ypredict == mydata[,phenoname])["TRUE"] / sum(table(ypredict == mydata[,phenoname]))
# print(paste("predictive power, proportion of right predictions=",myr2))
# }else{
# print ("assumed regression, continuous variable, change type flag to type classification if your response variable is discrete")
# myr2<-summary(lm(ypredict~mydata[,phenoname]))$r.squared
# print(paste("predictive power, R2=",myr2))
# }
#
#
# if(sink==T){
# sortedimportance<-data.frame(importance(rfall) )
# sortedimportance$r2<-myr2
# sortedimportance$var<-row.names(sortedimportance)
# # toreturn<-list(rfobject=rfall,importance=sortedimportance,myr2)
# write.tsv(sortedimportance,file = paste0("tables/",phenoname,"_",myfilename,"randomforest.tsv"))
# # print(sortedimportance)
# }
#
# return(rfall)
# }
