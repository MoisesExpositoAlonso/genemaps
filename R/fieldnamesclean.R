

cleanfieldphenos<-function(x){

gsub(paste(paste0("_",fieldcodes()), collapse = '|' ),
     "",
     x)

}

nonphenocolumns<-function(){
  c('id','name','country','latitude','longitude','kgroup')
}
