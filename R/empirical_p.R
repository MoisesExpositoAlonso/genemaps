
empirical_p<-function(value,distribution){
    # sorted<-sort(distribution)
    pemp<-as.numeric( table(distribution<value)['TRUE'] / length(distribution) )
  if(is.na(pemp)){pemp=0}
  return(pemp)
}
