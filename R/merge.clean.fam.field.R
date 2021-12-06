
merge.clean.fam.field<-function(fielddata,fam){
  fam$row=1:nrow(fam)
  fam[,'sample.ID']=fn(fam[,'sample.ID'])
  fielddata$id=fn(fielddata$id)

  Y=merge(
    fam,
    # fam[,'sample.ID'],
    by.x='sample.ID',
    fielddata,
    by.y='id',
    all.x=TRUE
  )
  Y$Fitness[is.na(Y$Fitness)]<-0

  return(Y)
}
