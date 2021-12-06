
fieldlocations<-function(site=NULL){
	fl<-data.frame(
			t=c(48.544886, 9.043042),
			m=c(40.408049, -3.835350))
	rownames(fl)<-c("lat","lon")
	if(is.null(site)){
		return(fl)
	}else if(site=='m'){
		return(fl$m)
	}else if(site=='t'){
		return(fl$t)
	}else{
		message("Provide one of the next c(m,t, or NULL)")
		return(fl)
	}
	
}