rastergrowingseason<-function(rasclim){

  tmpprec<-raster::subset(rasclim,paste0('prec',1:12))
  tmptmin<-raster::subset(rasclim,paste0('tmin',1:12))
  tmptmax<-raster::subset(rasclim,paste0('tmax',1:12))

  tmptmean<-((tmptmin + tmptmax)/2)/10
  # tmptmean<-((tmptmin + tmptmax)/2)

  gs <- (tmptmean >4) & (tmpprec > 2*tmptmean)

  names(gs)<-paste0('gs',1:12)

  gssum<-sum(gs)
  names(gssum)<-'gs_sum'


  for(i in 1:12){
    tmpprec[[i]][ gs[[i]]==0 ] <-NA
    tmptmin[[i]][ gs[[i]]==0 ] <-NA
    tmptmax[[i]][ gs[[i]]==0 ] <-NA
  }

  gs_prec= sum(tmpprec,na.rm=T)
  names(gs_prec)<-'gs_prec'
  # gs_cvprec= var(tmpprec,na.rm=T) / mean(tmpprec,na.rm=T)
  gs_tmin= min(tmptmin,na.rm=T)
  names(gs_tmin)<-'gs_tmin'
  gs_tmax=min(tmptmax,na.rm=T)
  names(gs_tmax)<-'gs_tmax'

  gsinfo<-stack(gs,gssum,gs_prec, gs_tmax,gs_tmin)
  return(gsinfo)
}
