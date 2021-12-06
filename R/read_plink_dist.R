
read_plink_dist<-function(file='../gemma/515g.dist.txt'){

  # read
  message("reading file ",file)
  d<-read.table(file,fill=T,header=T)

  # curate
  length(rownames(d))
  length(colnames(d))
  firstrow=rep(NA, ncol(d))

  d<-rbind(firstrow,d)

  rownames(d) <- colnames(d)[-ncol(d)] # the first and last no
  d<-d[,-ncol(d)]


  # fix columns
  colnames(d)<-rownames(d)<-
  gsub(colnames(d),pattern = "X", replacement = "", fixed = TRUE)

  # format nicer
  diag(d)<-0
  d[upper.tri(d)] <-d[lower.tri(d)]


  message('see extract of first 10 col/rows')
  print(d[1:10,1:10])


  return(d)
}


run_plink_dist<-function(){

  # plinkfile="../gemma/515g"
  #
  # # plink  --memory 32000 --bfile $plinkfile --distance --out $plinkfile
  #
  # dist=$plinkfile".dist"
  # iddist=$plinkfile".dist.id"
  # out=$plinkfile"dist.txt"
  #
  # awk '{print $1}' $iddist | tr "\n" "\t" >  $out
  # echo "\n" >> $plinkfile".dist"
  # cat  $dist >> $out

}
