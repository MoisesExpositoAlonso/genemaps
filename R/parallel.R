
startcluster<-function(no_cores=NULL){
require("parallel")
message("starting cluster for parallel operations")

# Calculate the number of cores
if(is.null(no_cores)){
  no_cores <- detectCores() - 1
  if(is.na(no_cores)) no_cores <- 2
  if(no_cores>4){no_cores=10}
}

 # Initiate cluster
# cl <- makeCluster(no_cores,type="FORK")
# Now we just call the parallel version of lapply, parLapply:
cl <- makeCluster(no_cores)

return(cl)
}

stopcluster<-function(cl){
  stopCluster(cl)
  message("cluster closed")

}
