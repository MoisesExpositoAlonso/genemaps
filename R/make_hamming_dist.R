
make_hamming_distances<-function(plinkfile='515g', path='../plink'){

#### for the raw count of number of SNPs different
# https://www.cog-genomics.org/plink2/distance#make_grm
# echo  "running plink to compute raw Hamming distance ..."
# plink  --memory 32000 --bfile $plinkfile --distance --out $plinkfile
# read_dist.R
fullname=paste0(path,'/',plinkfile)

message( "running plink to compute raw Hamming distance ...")
system(paste('plink  ',
             # '--memory 32000 ',
             '--bfile ',
             fullname,
             '--distance --out ',
             fullname)
       )

system(
       paste("awk '{print $1}'",paste0(fullname,".dist.id"), "| tr '\n' '\t' > ", paste0(fullname,".dist.tsv") )
       )
system(paste0("echo '\n' >> " ,fullname,".dist.tsv") )
system(paste0("echo '0' >> " ,fullname,".dist.tsv") )
system(paste0("cat ",fullname,".dist >> " ,fullname,".dist.tsv"))
message( "reading file " , paste0(fullname,".dist.tsv"), " ...")
# D<-read_relmat(paste0(fullname,".dist.tsv"))

D<-as.matrix(read.table(paste0(fullname,".dist.tsv"),fill=T,header=T))
dim(D)
rownames(D)=colnames(D)
stopifnot(nrow(D)==ncol(D))

diag(D) <- 0

return(D)
}
