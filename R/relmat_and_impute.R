read_relmat<-function(myfile){
  rawname<-tools::file_path_sans_ext(myfile)
  namerda<-paste0(rawname,".rda")

  if(file.exists(namerda)){
    load(namerda)
  }else{

  A<-as.matrix(read.table(myfile,fill=T,header=T))
  # A<-A[,-ncol(A)]
  rownames(A)=colnames(A)
  stopifnot(nrow(A)==ncol(A))
  save(A,file=namerda)

  }
  message("relationship matrix read")
  return(A)
}

make_relmat<-function(path="../reggenome", name="1001gregmapsmall.imputed"){
# #### Plink version (Relationship matrix rangin 0 to  2)                                                                              │1001_i09g01_maf0019relmat_2R.txt        1001_i09g01_maf0019relmat_2R.txt.RData  read_relmat.R
# echo "running plink to compute relationship matrix ..."                                                                              │> load("*relmat*")
# plink  --memory 32000 --file $plinkfile --make-rel --out $plinkfile                                                                  │1001_i09g01_maf0019relmat_2R.txt        1001_i09g01_maf0019relmat_2R.txt.RData  read_relmat.R
# echo "the relationship matrix file is:", $plinkfile".rel"                                                                            │> load("1001_i09g01_maf0019relmat_2R.txt.RData")
# echo "the column names of the relationship matrix file is:", $plinkfile".rel.id"                                                     │> ls()
#                                                                                                                                      │[1] "A" "K"
# # to make something readable in R                                                                                                    │> A[1:5,1:5]
# awk '{print $1}' $plinkfile".rel.id" | tr "\n" "\t" > $plinkfile"relmat_2R.txt"                                                      │           X88      X108      X139       X159    X265
# echo "\n" >> $plinkfile"relmat_2R.txt"                                                                                               │X88  1.5398000        NA        NA         NA      NA
# cat $plinkfile".rel" >> $plinkfile"relmat_2R.txt"                                                                                    │X108 0.0879091 1.7973600        NA         NA      NA
#                                                                                                                                      │X139 0.0853532 1.6695400 1.7299000         NA      NA
# Rscript read_relmat.R $plinkfile"relmat_2R.txt" # this generates an RData object with the matrix ready to use. Then in another script│X159 0.0674566 0.0252254 0.0255885 1.58639000      NA
#  you just have to load an object.

  fullname<-paste(sep="/", path, name)
    if(file.exists(paste0(fullname,".relmat.tsv"))){
    A<-read_relmat(paste0(fullname,".relmat.tsv"))
  }else{

  system(paste("plink","--file", fullname, "--make-rel --out", fullname ))
  system(
         paste("awk '{print $1}'",paste0(fullname,".rel.id"), "| tr '\n' '\t' > ", paste0(fullname,".relmat.tsv") )
         )
  system(paste0("echo '\n' >> " ,fullname,".relmat.tsv") )
  system(paste0("cat ",fullname,".rel >> " ,fullname,".relmat.tsv"))
  A<-read_relmat(paste0(fullname,".relmat.tsv"))
  }

  return(A)
}


#####**********************************************************************#####
read_kinship<-function(myfile){

  rawname<-tools::file_path_sans_ext(tools::file_path_sans_ext(myfile))
  namerda<-paste0(rawname,"kinship.rda")

   if(file.exists(namerda)){
    load(namerda)
  }else{
  K<-as.matrix(read.table(myfile,fill=T,header=F))

  fam=read.table(paste0(rawname,".fam"))
  colnames(K) <-rownames(K) <- moiR::fc(fam[,1])

  save(file=paste0(rawname,"kinship.rda"),K)

  }
  message("kinship matrix read")
  return(K)
}

make_kinship<-function(path="../reggenome", name="1001gregmapsmall.imputed"){
# ##### Emmax version (Kinship matrix ranging 0 to 1. to transform to Relationship matrix multiply by two)
# outplink="1001gtoemmax"
#
# echo "running emmax to compute the kinship matrix ..."
# plink --memory 32000 -bfile $plinkfile --output-missing-genotype 0 --recode transpose 12  --out $outplink
# emmax-kin -v -h -s -d 10 $outplink
# echo "the kinship matrix file is:", $outplink".hIBS.kinf"
#
# Rscript read_kinship.R $outplink".hIBS.kinf" # this is going to produce an RData object with the kinship matrix. In this case we do not have the rownames but they are the same order as in the plink .fam file

  fullname<-paste(sep="/", path, name)
  if(file.exists(paste0(fullname,".hIBS.kinf"))){
    K<-read_kinship(paste0(fullname,".hIBS.kinf"))
  }else{

  system(paste("plink","--file", fullname, "--output-missing-genotype 0 --recode transpose 12  --out", fullname ))
  system(paste("emmax-kin -v -h -s -d 10 ", fullname))
  # system(
  #        paste("awk '{print $1}'",paste0(fullname,".rel.id"), "| tr '\n' '\t' > ", paste0(fullname,".relmat.tsv") )
  #        )
  # system(paste0("echo '\n' >> " ,fullname,".relmat.tsv") )
  # system(paste0("cat ",fullname,".rel >> " ,fullname,".relmat.tsv"))

    K<-read_kinship(paste0(fullname,".hIBS.kinf"))
  }

  return(K)
}




#####**********************************************************************#####
make_Ainv<-function(A, mirrorupper=T, cleannameX=T){
  if(mirrorupper){
    A[upper.tri(A)]<-A[lower.tri(A)]
  }
  genonames<-colnames(A)
  if(cleannameX){   gsub("X","",genonames)  }

  A<-as.matrix(A)

  Ainv<-as(MASS::ginv(A), "dgCMatrix")
  rownames(Ainv) = rownames(Ainv) = colnames(A)
  return(Ainv)
}

