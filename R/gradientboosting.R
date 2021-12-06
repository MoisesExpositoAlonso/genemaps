#' Gradient boosting wrapper
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


def.fitcontrol<-function(){
  fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)
  return(fitControl)
}

grandientboosting<-function(myformula, dat, tuningmode=FALSE,...){
set.seed(1)

# From Caret package: https://topepo.github.io/caret/model-training-and-tuning.html
# number of iterations, i.e. trees, (called n.trees in the gbm function)
# complexity of the tree, called interaction.depth
# learning rate: how quickly the algorithm adapts, called shrinkage
# the minimum number of training set samples in a node to commence splitting (n.minobsinnode)



fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)

if(tuningmode==TRUE){

modtune <-  expand.grid(interaction.depth = c(1, 5, 10),
                        n.trees = (1:30)*50,
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
modtune <- train(myformula, data = dat,
                 method = "gbm",
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 ## Now specify the exact models
                 ## to evaluate:
                 tuneGrid = gbmGrid,
                 ...)
trellis.par.set(caretTheme())

ptune<-plot_grid(
  ggplot(gbmFit2metric = "Rsquared"),
  ggplot(gbmFit2,metric='RMSE')
)

bestmodel<-best(gbmFit2$results,metric='RMSE',maximize = FALSE)

toreturn<-list(bestmodel, mod=modtune,plot=ptune)

}else{
mod <- train(myformula, data = dat,
                 method = "gbm",
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 ...)

toreturn<-list(mod)

}

return(toreturn)

}




window_gb<-function(){

}


#
#
# mydata
#
# grandientboosting<-function(myformula,dat, tuningmode=TRUE){
#
# # From Caret package: https://topepo.github.io/caret/model-training-and-tuning.html
# # number of iterations, i.e. trees, (called n.trees in the gbm function)
# # complexity of the tree, called interaction.depth
# # learning rate: how quickly the algorithm adapts, called shrinkage
# # the minimum number of training set samples in a node to commence splitting (n.minobsinnode)
#
#
# gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 10),
#                         n.trees = (1:30)*50,
#                         shrinkage = 0.1,
#                         n.minobsinnode = 20)
# }
#
# nrow(gbmGrid)
#
# set.seed(1)
# gbmFit2 <- train(myformula, data = training,
#                  method = "gbm",
#                  trControl = fitControl,
#                  ## This last option is actually one
#                  ## for gbm() that passes through
#                  verbose = FALSE,
#                  ## Now specify the exact models
#                  ## to evaluate:
#                  tuneGrid = gbmGrid)
#
# bestmodel<-best(gbmFit2$results,metric='RMSE',maximize = FALSE)
#
# gbmFit2$results[bestmodel,1:6]
#
# whichTwoPct <- tolerance(gbmFit2$results, metric = "best",maximize = TRUE)
#
#
# gbmFit3$results[whichTwoPct,1:6]
#
#
# trainControl
#
# trellis.par.set(caretTheme())
# plot(gbmFit2)
# plot(gbmFit2, metric = "Rsquared")
#
#
#
# library(caret)
# fitControl <- trainControl(## 10-fold CV
#                            method = "repeatedcv",
#                            number = 10,
#                            ## repeated ten times
#                            repeats = 10)
#
#
# gbmFit1 <- train(myformula, data = training,
#                  method = "gbm",
#                  trControl = fitControl,
#                  ## This last option is actually one
#                  ## for gbm() that passes through
#                  verbose = FALSE)
# gbmFit1
#
#
#
# gradientboostingcv<-function(){
#
#
# }
#
#
# randomForestcv<-function(mydata,var,mypredictors,cross=5,seed=1, type="regression",sink=T, out='tmptables' ,filename="",verbose=T){
#
# # checking arguments
# if(!type %in%  c('regression','classification')){message('provide an appropriate random forest type')}
#
# # load llibraries
# library(caret)
# library(randomForest)
#
# # generate the cross validation partitions
# if(verbose) message(paste('creating',cross, 'partitions of the data ...'))
# group<-createFolds(mydata[,'beta'] , k=cross,list=FALSE)
#
# # set up several parameters
# allpredicted=c()
# varcol=which(colanmes(mydata)==var)
#
# if(verbose) message(paste('running Random Forest...'))
# rfmodels<-lapply(1:cross,FUN =
#           function(l){
#             message(paste('# ',l))
#
#             train <- mydata[group != l,]
#             test <- mydata[group == l,-varcol]
#
#             if(type=="classification"){ mydata[,var] = factor(mydata[,var])
#             }else{mydata[,var] = as.numeric(mydata[,var])}
#
#             myformula<-formula(paste(var," ~ ", paste(mypredictors,collapse = '+')))
#
#             fit <- train(myformula, data=mydata,method='rf', prox=TRUE )
#             # fit <- randomForest(myformula, data=mydata,ntree=500,na.action=na.omit,keep.inbag = TRUE,importance=TRUE) # removed proximity
#
#             mypred <-predict(fit,test)
#             allpredicted<-c(allpredicted,mypred)
#
#           return(fit)
# } )
#
# # combine all random forest
# if(verbose) message("merging the random forest")
# rfall<-randomForest::combine(rfmodels[[1]],rfmodels[[2]],rfmodels[[3]],rfmodels[[4]],rfmodels[[5]])
#
# # Generate accuracy estimation
# if(type=="classification"){
#     myr2= table(ypredict == mydata[,phenoname])["TRUE"] / sum(table(ypredict == mydata[,phenoname]))
#     if(verbose) message(paste("predictive power, proportion of right predictions=",myr2))
# }else{
#     myr2<-summary(lm(allpredicted~mydata[,phenoname]))$r.squared
#     if(verbose) message(paste("predictive power, R2=",myr2))
# }
#
# # sink info
# if(sink==T){
#     sortedimportance<-data.frame(importance(rfall) )
#     sortedimportance$r2<-myr2
#     sortedimportance$var<-row.names(sortedimportance)
#     write.tsv(sortedimportance,file = paste0(out, '/',var,"_",filename,"randomforest.tsv"))
# }
#
# return(rfall)
# }
#
