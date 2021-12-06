get_realfit<-function(fitness='Fitness_mlp',relative=TRUE){
  Y=merge(genomes$fam,by.x="sample.ID", dry[,c("id",fitness)], by.y="id",all.x=T)[,7]

  Y[is.na(Y)]<-mean(Y,na.rm=TRUE)

  if(relative) {
    Y<-relative(Y)
  }
  return(Y)
}

get_topSNPs<-function(p=100,navalue=0.5){

    ## Best positions
    data('gwamlp')
    gwares<-gwamlp

    maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
    maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

    subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[p+1]),]
    subg$rs<-paste(subg$chr,subg$pos,sep='_')
    head(subg)

    ###############################################################################
    ## Subset of matrix
    data("genomes")
    G <- genomes$genotypes
    Map <- genomes$map
    dim(Map)
    head(Map)
    Map$rs<-paste(Map$chr,Map$physical.pos,sep='_')
    Map$order<-1:nrow(Map)

    matched<-merge(Map, subg, by='rs')
    dim(matched)
    head(matched)

    X=change01(inputena.mat(G[,matched$order],value = navalue))

    return(X)
}

top_cols<-function(genomes,p=100,gwa='gwamlp'){

  ## Best positions
  eval(parse(text = paste0("data(",gwa,")")))
  eval(parse(text = paste0("assign('gwares', ",gwa,")")))

  maxgwa<-dplyr::filter(gwares,p_lrt == min(p_lrt))
  maxgwa<-dplyr::filter(gwares,abs(beta) == max(abs(beta)))

  subg<-gwares[which(gwares$p_lrt < sort(gwares$p_lrt)[p+1]),]
  subg$rs<-paste(subg$chr,subg$pos,sep='_')
  head(subg)

  ###############################################################################
  ## Subset of matrix
  G <- genomes$genotypes
  Map <- genomes$map
  dim(Map)
  head(Map)
  Map$rs<-paste(Map$chr,Map$physical.pos,sep='_')
  Map$order<-1:nrow(Map)

  matched<-merge(Map, subg, by='rs')
  dim(matched)
  head(matched)

  return(matched$order)
}

get_genomematrix<-function(size=1000,positions=NULL,type=c('random','window','top'),start=1){

# data(genomes)
G <- genomes$genotypes

if(!is.null(positions) ){
  message('subsetting the provided positions')
  stopifnot(is.numeric(positions))
  X = inputena.mat(G[,positions])
}else if(type=='random'){
  message('subsetting ', size, ' random positions ')
  X = inputena.mat(G[,sample(1:ncol(G),size = size)])
}else if(type=='window'){
  message('subsetting a window of ', size, ' base pairs starting at ',start)
  X = inputena.mat(G[,start:(start+size)])
}else if(type=='top'){
  X<-get_topSNPs()
}else{
  stop('None of the possibilities were given')
}

return(X)
}


# X=matrix(ncol=5,nrow=10)
# X[is.na(X)]<-sample(size=length(X),x = c(0,1),replace = TRUE)
# sparsity=0

simupheno<-function(X,meanN=1,sdN=5,sparsity=0,type="additive", epistasis=0.5,heritability=1, zerobound=TRUE){

  # Simulate Gaussian effects
  eff<-rnorm(n=ncol(X),meanN,sdN)

  # Apply sparsity
  eff<- eff * rbinom(n=ncol(X),size = 1,prob = 1-sparsity)

  if(type=='additive'){
    # Get the phenotype based on effect of SNPs
    Y= X %*% eff

  }else if(type=='multiplicative'){
    Y= apply(X,1, function(i) prod(1+ i* eff))

  }else if(type=="epistatic"){
    # if epistatic add multiplicative terms with an extra term
    effx= tcrossprod(eff, eff)/2
    effx = effx * epistasis
    extraeffect= apply(X,1,function(i) sum( tcrossprod(i) * effx ))

    Y=Y+extraeffect
  }else{
    Stop('Provide type additive, multiplicative or epistatic!')
  }

  # Reducte heretability from 1
  if(heritability!=1){
    Y=addnoise(Y,noise = 1-heritability)
  }

  if(zerobound==TRUE){
    Y[Y<0] <-0
  }

  return(list(Y=Y,eff=eff))
}


reduceH2 <- function(Y,rho) {
  # from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
n     <- length(Y)                    # length of vector
rho   #<- 0.6                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- Y        # fixed given data
x2    <- rnorm(n, 2, 0.5)      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector

cor(x1, x)
return(x)
}


buildpheno<-function(X,s,type='multiplicative')
{
  if(type=='multiplicative'){
    y=apply(X,1,function(x) prod(1+x*s))
  }else if(type=='additive'){
    y=apply(X,1,function(x) 1+sum(x*s))
  }else{
    stop('types implemented are multiplicative and additive')
  }
  return(y)
}
