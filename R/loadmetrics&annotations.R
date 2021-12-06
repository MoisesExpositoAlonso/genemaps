load_metrics_n_annotations<-function(){

  load("dataint/gmetrics.rda")
  metricscols<-c("maf", "fst" ,"pi" ,"taj" ,"sweep_lr", "sweep_alpha")

  ### Load data - an
  load('dataint/an01.rda')
  ancols<-c("intergenic", "intron", "UTR3", "UTR5" ,"exon", "synonymous", "nonsynonymous" ,"exon_noncoding")

  allmetrics<- cbind(gmetrics, an)
  rownames(allmetrics)<-allmetrics$SNP
  allmetcols<-c(metricscols,ancols)



  assign("gmetrics", gmetrics, envir=globalenv())
  assign("metricscols", metricscols, envir=globalenv())
  assign("an", an, envir=globalenv())
  assign("ancols", ancols, envir=globalenv())
  assign("allmetrics", allmetrics, envir=globalenv())
  assign("allmetcols", allmetcols, envir=globalenv())


}
