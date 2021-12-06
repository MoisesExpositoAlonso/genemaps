rasterevapo<-function(r,climname){
  require(raster)

filename<-paste0('dataint/',climname,'evapo.grd')

# if(file.exists(filename)){
#   message('evapotranspiration raster file found for this climate layer')
#   et <- stack(filename)
#
# }else{

  # https://duncanjg.wordpress.com/2012/11/29/building-more-informative-climate-layers-for-species-distribution-modelling/
  # potential evapotranspiration

  require("dismo")
  require("EcoHydRology")

  names(r)<-c(paste("mint",1:12,sep=""),paste("maxt",1:12,sep=""),paste("prec",1:12,sep=""))

  lat_rad <- coordinates(r)[, 2] * pi/180
  length(lat_rad)

  for (i in 1:12) {
    evap <- raster(r, 1)
    Tmax <- values(subset(r, i + 12))/10
    Tmin <- values(subset(r, i))/10
    d <- data.frame(day = (30 * i) - 15, Tmin, Tmax, lat_rad)
    d[is.na(d)] <- 0
    Es_PET <- PET_fromTemp(Jday = d$day, Tmax_C = d$Tmax, Tmin_C = d$Tmin, lat_radians = d$lat_rad) * 1000
    values(evap) <- Es_PET
    if (i == 1) {
        PET <<- brick(evap)
    }
    if (i > 1) {
        PET <<- addLayer(PET, evap)
    }
  }

  months<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  names(PET)<-months


  Bucket<-raster(PET,1)

  for (n in 1:2) {
    for (i in 1:359) {
        mn <- 1 + i%/%30
        NewAET <- raster(PET, 1)
        NewBucket <- values(Bucket)
        rain <- values(raster::subset(r, 24 + mn))/30
        alpha <- (NewBucket - 200)/300
        evap <- values(raster::subset(PET, mn)) * alpha * 0.8  ##A fudge factor for stomatal control.
        NewBucket <- NewBucket + (rain) - evap
        NewBucket[NewBucket > 500] <- 500
        NewBucket[NewBucket < 200] <- 200
        values(Bucket) <- NewBucket
        values(NewAET) <- evap * (NewBucket > 200)
        if (n > 1 && (i%%30) - 15 == 0) {
            if (mn == 1) {
                AET <<- brick(NewAET)
            }
            if (mn > 1) {
                AET <<- addLayer(AET, NewAET)
            }
        }
    }
  }

  names(AET)<-paste0('aet',1:12)
  names(PET)<-paste0('pet',1:12)
  et<-stack(AET,PET)

  writeRaster(file=filename,et)
# }

return(et)
}
