
findSNPname <- function(i,snptable,ann){

  print(i)
  annsub<-subset(ann,ann$V1 == paste("Chr",snptable[i,"chr"],sep='') & ann$V4<snptable[i,"pos"] & ann$V5> snptable[i,"pos"]  )
  annsub_select<-subset(annsub,annsub$V3 %in% c('gene'))
  annsub_select_an<-as.character.factor(annsub_select$V9)
  if(length(annsub_select_an) != 0L){
  genename<-strsplit(annsub_select_an,split = 'Name=')[[1]][2]
  # annsub_collapse<-paste(annsub$V9,collapse = "|")
  # annsub_collapse_1<-gsub(annsub_collapse,pattern = "Name=",replacement = "")
  # annsub_collapse_2<-gsub(annsub_collapse_1,pattern = "ID=",replacement = "")
  # annsub_collapse_3<-gsub(annsub_collapse_2,pattern = "Parent=",replacement = "")
  # strsplit(annsub_collapse_3,split = "|",fix=T)[[1]]
  # grep(pattern = "AT[1-5]G*",x = strsplit( annsub_collapse,split = c(";","=","|") ) )
  # genename<-as.character.factor(annsub$V9[grep("ID=",annsub$V9)] )
  # genename_list<-strsplit(genename,split = ";",fixed = T)[[1]]
  # genename_list_parse<-genename_list[grep("ID=",genename_list)]
  # genename_list_parse_clean<-sub(x = genename_list_parse,pattern = "ID=",replacement = "")
  # newsnptable$genenames[i]<-genename_list_parse
  return(genename)
  }else{
    print('not found')
    return(NA)
    }

}

samplebackground<-function(n,ebioroot='~/ebio/'){
  
  message("reading the TAIR annotations file")
  ann<-read.table(paste(ebioroot,'abt6_projects8/ath_1001G_history/TAIR10_GFF3_genes_transposons.gff.txt',sep=''))
  ann<-subset(ann,V3=='gene')
  ann<-as.character.factor(ann$V9)
  genename<-sapply(ann, function(x) strsplit(x,split = 'Name=')[[1]][2])
  genename<-moiR::fc(genename)
    
  return(sample(genename,n))
}


findgenename<-function(snptable=genelist,ebioroot='~/ebio/') {
  
  if(colnames(snptable)[1] != 'chr' | colnames(snptable)[2] != 'pos')  {
    stop(' please rename your snptable as two columns: chr pos')
        }
  
  message("reading the TAIR annotations file")
  ann<-read.table(paste(ebioroot,'abt6_projects8/ath_1001G_history/TAIR10_GFF3_genes_transposons.gff.txt',sep=''))
  head(ann)
  
  message("looking for genes hit by the SNPs")
  genenames=sapply(1:nrow(snptable), function(i) findSNPname(i,snptable,ann) )
  
  snptable$genenames=genenames
  return(snptable)
}

findgenename_SNPeff<-function(snptable=genelist,ebioroot='~/ebio/') {
  
  #### check for the column names
  if(colnames(snptable)[1] != 'chr' | colnames(snptable)[2] != 'pos')  {
    stop(' please rename your snptable as two columns: chr pos')
        }
  
  #### write tmp file to call the python script
  write.table(file = ".tmp.tsv", sep="\t" , snptable,col.names = F, row.names=F)
  tmpfile= paste0(getwd(), "/.tmp.tsv")
  
  command <- paste0("python ",ebioroot,"abt6_projects9/ath_1001g_foreverybody/minivcftools/findSNPeff_fromlist.py ", tmpfile)
  system(command)
  
  anfile= paste0(".tmp.tsv_genearound")
  anadd=read.table(anfile,header=F)
  return(anadd)
}

.parser_SNPeff <-  function(x){
tmp<- strsplit( strsplit(as.character(x),split = "EFF=")[[1]][2] ,split =  "(" ,fixed = T) [[1]] [1]
return(as.character(tmp))
}

parser_SNPeff <- function(anadd){
  as.character(sapply(anadd[,3], FUN = .parser_SNPeff) )
}

genes_SNPeff <- function(anadd){
  tmp<- sapply(anadd[,3],.genes_SNPeff )
  # print(tmp)
  # tmp2 <- paste0(as.character(tmp), collapse = " | ")
  # print(tmp2)
  # return(tmp2)
  # reduce(tmp,paste(collapse = "|"))
  # return(as.character(tmp))
  sapply(tmp,FUN = function(x) paste(x, collapse = " "))
}
.genes_SNPeff <- function(x){
x=as.character(x)
splitted<- strsplit(x, split = "|",fixed=T)
locations<- grep(pattern = "AT[1-5]G",perl = T, splitted[[1]])
parsed<- splitted[[1]] [locations]
parsed2 <-gsub(parsed,pattern = "\\..",fixed = F, replacement = "") # this replaces a dot (\\.) and any character afterwards(.). The backslash is necessary so it knows it is an actual dot.
parsed2 <- unlist(unique(parsed2))
return(parsed2)
}

loadgoslim<-function(){
# gocat<-read.table('~/ebio/abt6_projects8/ath_1001G_history/GO_Ath/TAIR_GO_slim_categories.txt',fill=T,header=T)
# head(gocat)
print("loading goslim datasets stored at:")
filegoslim<-'~/ebio/abt6_projects8/ath_1001G_history/GO_Ath/ATH_GO_GOSLIM.txt'
print(filegoslim)
goslim<-read.delim(filegoslim,fill=T,header=F)
head(goslim)


goslim$V1 # at gene name
goslim$V2 # locus id
goslim$V3 # at gene name or splicing form
levels(goslim$V4) # preambule
levels(goslim$V5) # word description
levels(goslim$V6) # actual go term
(goslim$V7) # tair keyword
(goslim$V8) #  F=molecular function, C=cellular component, P=biological process
levels(goslim$V9)  # classification description
# (goslim$V10)

# source("https://bioconductor.org/biocLite.R")
# biocLite("GEOsearch")

# https://www.arabidopsis.org/servlets/Search?action=new_search&type=keyword
# https://www.arabidopsis.org/servlets/Search?action=new_search&type=gene
# https://www.arabidopsis.org/servlets/Search?type=general&action=new_search
# return(goslim)
assign("goslim", goslim, envir=globalenv())
}

checkgoslim <- function(){
  if( !("goslim" %in% ls(envir=globalenv()))) {
  message("goslim not found, loading...")
  loadgoslim() 
}else{
  # message("goslim found")
  } 
}

humanname <- function(genes,columngenes="genenames", goodcolumns=c("genenames","chr","pos","V5","V9")){
checkgoslim()

if(is.null(columngenes)){message("you need to provide a data frame for genes, and the column name that you want to merge with the gene names, normally V5 and V9 that are annotations")}

#### merge with goslim annotations
merged<-merge(genes,by.x=columngenes,goslim,by.y="V1")

#### get important column and remove duplicates
namesan.unique=unique(merged[,goodcolumns])

#### summarize the 

go_broad <- as.matrix(tapply(namesan.unique$V9,namesan.unique$genenames,FUN = function(x){paste(unique(x),collapse = " | ")}))
go_broad <- data.frame(genenames=rownames(go_broad), go_broad=go_broad)
go_narrow <-as.matrix(tapply(namesan.unique$V5,namesan.unique$genenames,FUN = function(x){paste(unique(x),collapse = " | ")})) 
go_narrow<- data.frame(genenames=rownames(go_narrow), go_narrow=go_narrow)

newgenelist_go <-left_join(newgenelist, go_broad,by=c("genenames")) %>% left_join (.,go_narrow,by=c("genenames"))

print("the resulting table has the next structure")
print(str(newgenelist_go))

return(newgenelist_go)

}

humanname_specific <- function(genes,columngenes="genenames"){

library(org.At.tair.db)
library(DBI)

geneinfolist<-select(org.At.tair.db, keys=na.omit(genes[,columngenes]),columns=c("SYMBOL","GENENAME"))
geneinfolist <- unique(geneinfolist)


geneinfolist_tmp<-as.matrix(tapply(geneinfolist$SYMBOL,geneinfolist$TAIR,FUN = function(x){paste(unique(x),collapse = " | ")}))
geneinfolist_wordy_tmp<-as.matrix(tapply(geneinfolist$GENENAME,geneinfolist$TAIR,FUN = function(x){paste(unique(x),collapse = " | ")}))

newgeneinfolist<-data.frame(TAIR=row.names(geneinfolist_tmp),SYMBOL=geneinfolist_tmp[,1],DESCRIP=geneinfolist_wordy_tmp)

print("the resulting table has the next structure")
print(str(newgeneinfolist))

return(newgeneinfolist)

}

grepgenes<-function(genes=myhits,togrep="port",columngo="V9",type="union",wannaprint=F,columngenes="genenames"){
checkgoslim()
  
merged<-merge(genes,by.x=columngenes,goslim,by.y="V1")
head(merged)
dim(merged)

if(type=="union"){
found<-unique(unlist(sapply(togrep,FUN=function(x) grep(x, merged$V5))))

}else if(type=="intersection"){
found1<-unlist(sapply(togrep,FUN=function(x) grep(x, merged$V5)))
found<-found1[duplicated(found1)]

}else{ print ("need to provide intersection or union in type flag!")}

# message("returned all SNPs merged and those grepped by the keyword")
found<-merged[found, ] 
return(list(merged=merged,found=found))
# return(searched)
}

count_grep<-function(genes=myhits,togrep="port",columngo="V9",type="union",wannaprint=F, columngenes="genenames",searched=NULL){

if(is.null(searched)){
searched=grepgenes(genes=genes,togrep=togrep,columngo=columngo,type=type,wannaprint=wannaprint,columngenes=columngenes)
# searched<-merged$V5[found ] 
}

annot=searched$found[,columngo]


total_genes<-length(na.omit(unique(genes[,columngenes])))
total_merged <- nrow(searched$merged)
total_unique_merged <- length(unique(searched$merged[,1]))

wordhits<-length(annot)
uniquegeneword <- length(unique(searched$found[,1]))

results<-list( unique_word_genes= uniquegeneword,
                    all_word_found=wordhits,
                    unique_merges= total_unique_merged,
                    all_merges=total_merged,
                    all_genes_provided= total_genes
                    )

if(wannaprint==T){print(results) }

return(results)
}


plotcounts<-function(test,neutral,name,...){
require(cowplot)
require(ggplot2)
p<-ggplot(data.frame(hits=neutral))+geom_histogram(aes(x = hits),fill="black",color="white",...) + ggtitle(name)+ geom_vline( aes(xintercept= test),color="red" )
print(p)
}

genepoisson <- function(test,neutral){
  
print(paste("average # found in neutral sets: ",mean(neutral)) )
print(paste("total # genes in my hits: ",test) )
print (paste("the p-value under Poisson:", ppois(test,lambda = mean(neutral),lower.tail = F) ) )

}

genetstudent <- function(test,totaltest,neutral,totalneutral){
  
print(paste("average proportion in neutral sets: ", mean(neutral/totalneutral)) )
print(paste("total # genes in my hits: ",mean(test/totaltest)) )
print (paste("the p-value under t distribution:", t.test(neutral/totalneutral, mu= test/totaltest,alternative="less" )$p.val ) )

}


neutral_go<-function(times=10,togrep="port",size=30,sizem=size*2,category=NULL,type="union",columngo ="V9"){
print(paste("generating empirical distribution of",times, "random draws...") )
  
  thecounts<- sapply(1:times, FUN = function(x){
    
    # message(dim(goslim))
    # message(size)
    
    #### first sample a bit higher in case some are not found
    # samplenames<-sample(as.character.factor(goslim[,'V1']), size=sizem,replace = F)
    samplenames<-sample(as.character.factor(goslim[,'V1']), size=size,replace = F)
    
    # searched <- grepgenes(genes=data.frame(genenames=samplenames),togrep = togrep,type = type,columngo = columngo )
    # searched$merged$genenames <- as.character(searched$merged$genenames)
    # searched$found$genenames <- as.character(searched$found$genenames)
    # 
    # unique(searched$merged$genenames)
    # unique(searched$found$genenames)
    # 
    # #### secondly subsample. 
    # subsample<- sample(searched$merged$genenames, size = size)
    # 
    # searched$merged <-filter(searched$merged, genenames %in% subsample)
    # searched$found <- filter(searched$found, genenames %in% subsample)
    
      
    found<-count_grep(genes = data.frame(genenames=samplenames),wannaprint = F )#,searched = searched)
    
    if(!is.null(category)){counted<-found[[category]]}
    else{counted<-found}
    
    
    return(as.numeric(counted))
  }
  )
  
return(thecounts)
}



wrap_annotation_test<-function(genes,togrep,type="union",columngo="V5",wannaprint=T,times=1000){

message("searching for type: ",type)
message("of the keywords: ", paste(togrep,collapse= " "))

searched <- grepgenes(genes = genes,columngo = columngo,togrep = togrep,wannaprint = T,type = type)
# unique(as.matrix(searched$found$V5))
message(unique(as.matrix(searched$found$V5)))

message("THE FOCAL GENE SET COUNT")
counted<- count_grep(genes = genes,wannaprint = T,columngo = columngo,searched = searched,columngenes = "genenames")
coun=as.numeric(counted)

message("THE NEUTRAL GENE SET COUNT")
neut <- neutral_go(times = times,togrep = togrep,columngo = columngo,  size = counted$all_genes_provided,  type = type)
message(paste("average size of random sample:", mean(neut[3,])))

genepoisson(coun[1],neutral = neut[1,] )
genetstudent(coun[1],coun[3],neut[1,],neut[3,] )

return(searched)
}
