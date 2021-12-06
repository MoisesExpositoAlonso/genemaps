forceNC<-function(x,type){
	if(type %in% c("interger","numeric","matrix")){
		x<-fn(x)
	}else if(type %in% c("character")){
		x<-fc(x)
	}else{
		stop("Unknown type")
	}
}
	
cleandatasettypes<-function(x){
	
	apply(data.frame(x),2, function(i) forceNC(i, class(i)))
	
}