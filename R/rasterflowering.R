
rasterflowering<-function(rasclim,climname,force=F){
require(raster)

filename<-paste0('dataint/',climname,'lifespan.grd')

# if(file.exists(filename) & force==F){
#   message('lifespan file found for this climate layer')
#   lifespan <- stack(filename)
#
# }else{

  acc515env<-readRDS('dataint/acc515env.rda')[euroacc,]
  newbioclim<-stack('dataint/newclim.grd')
  # broadeuroextents<-broadeuroextents()
  # euroacc <- which( (( fn(acc515$latitude) > broadeuroextents$ylim[1] & fn(acc515$latitude) < broadeuroextents$ylim[2] ) &( fn(acc515$longitude) > broadeuroextents$xlim[1] & fn(acc515$longitude) < broadeuroextents$xlim[2] ) ) ==T)

  require(randomForest)

  ####  Extract info from field for mean flowering
  load('dataint/flowermean.rda')
  flowermean<-flowermean[euroacc,]

  flowermean<-flowermean[, which(colnames(flowermean) %in% c('meanFlowering', colnames(acc515env)) )]
  flowermean<-na.omit(flowermean)
  dim(flowermean)

  frf<-randomForest(y= flowermean$meanFlowering , x= flowermean[,colnames(acc515env)] )
  frf
  predictedflor<-raster::predict(model=frf,object=newbioclim)

  ####  Extract info from field for variance flowering
  load('dataint/flowervar.rda')

  flowervar<-flowervar[, which(colnames(flowervar) %in% c('varFlowering', colnames(acc515env)) )]
  flowervar<-na.omit(flowervar)
  dim(flowervar)

  frf<-randomForest(y= sqrt(flowervar$varFlowering) , x= flowervar[,colnames(acc515env)] )
  frf
  predictedvar<-raster::predict(model=frf,object=newbioclim)

  ####  put together twopredicted rasters
  lifespan<-stack(predictedflor, predictedvar)
  names(lifespan)<-c('cyclelength','cyclevar')


  # writeRaster(file=filename,lifespan)
# }

return(lifespan)
}


flowering_vs_growingseason<-function(rasclim){
  ####  generate derived layers
  cyclegs<-rasclim$gs_sum - (rasclim$cyclelength / 30)

  rasclim<-addLayer(rasclim,cyclegs)
  names(rasclim)[nlayers(rasclim)] <- 'cyclevsgs'
  return(rasclim)
}
