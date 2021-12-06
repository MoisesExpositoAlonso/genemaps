# extractclimate<-function(myrasterstack,coords){
#
# if(any(apply(coords, 2,class) != c('numeric','numeric') )){
#   message('Coords must be numeric!')
#   message('Attempting transformation...')
#   coords= apply(coords,2,as.numeric)
#   }
#
# message('Extracting values from raster')
# climgenom=
#   sapply(1:19, function(b){
#   raster::extract(myrasterstack[[b]],coords)
#
# })
# colnames(climgenom)<-names(myrasterstack)
#
# return(climgenom)
# }
#
#
#
# cropenvironment <- function(mybioclim,xlim=c(-10.5,+ 53),ylim=c(32,65),replace=F , addPCA=F)  {
#
#   Range=extent(c(xlim,ylim))
#   EuropeClim = crop(mybioclim, Range)
#
#   if(addPCA==T){
#     EuropeClim<-morelayers_2_bioclim(bioclim = EuropeClim,newlayerslist = stack(list(PC1=PC1,PC2=PC2,PC3=PC3)) )
#   }
#
# return(EuropeClim)
# }
