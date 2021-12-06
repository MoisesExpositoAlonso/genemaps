inputena.mat<-function (X,value=0.5){
  X[is.na(X)]<-value
  return(X)
}

meancent.mat<-function (X.)
{
  X = apply(X., 2, function(x) {
    mu = mean(x)
    x - mu
  })
  X
}

meanvarcent.<-function(x){
  x=fn(x)
  (x - mean(x, na.rm = T))/sd(x, na.rm = T)
}
meanvarcent<-function (x,row_or_col=2 ){

  if(class(x) == "numeric"){
    return(meanvarcent.(x))
  }else if(class(x)=='matrix'){
    return(apply(x,row_or_col, meanvarcent.))
  }else{
    stop('Data type not numeric nor matrix')
  }
}
meanvarcent.mat<-meanvarcent

narm.mat<-function (X)
{
  X = apply(X, 2, function(x) {
    x[x == (-1)] <- 1
    x
  })
  X = apply(X, 2, function(x) {
    x[is.na(x)] <- 1
    x
  })
  X
}
logistic<-function (x, L = 1, x0 = 0.05, k = 10)
{
    L/(1 + exp(-k * (x - x0)))
}

logit<-function (p, percents = range.p[2] > 1, adjust)
{
    range.p <- range(p, na.rm = TRUE)
    if (percents) {
        if (range.p[1] < 0 || range.p[1] > 100)
            stop("p must be in the range 0 to 100")
        p <- p/100
        range.p <- range.p/100
    }
    else if (range.p[1] < 0 || range.p[1] > 1)
        stop("p must be in the range 0 to 1")
    a <- if (missing(adjust)) {
        if (isTRUE(all.equal(range.p[1], 0)) || isTRUE(all.equal(range.p[2],
            1)))
            0.025
        else 0
    }
    else adjust
    if (missing(adjust) && a != 0)
        warning(paste("proportions remapped to (", a, ", ", 1 -
            a, ")", sep = ""))
    a <- 1 - 2 * a
    log((0.5 + a * (p - 0.5))/(1 - (0.5 + a * (p - 0.5))))
}

inv.logit<-function (f, a = 0.025)
{
    a <- (1 - 2 * a)
    (a * (1 + exp(f)) + (exp(f) - 1))/(2 * a * (1 + exp(f)))
}

logrelative<-function (f, min = 1)
{
    f <- log10(f + min)
    tmp <- f/mean(f, na.rm = T)
    return(tmp)
}

normalize<-function (x)
{
    if (class(x) != "numeric") {
        message("non-numeric data provided")
        message("attempting to force numeric")
        x = fn(x)
    }
    (x - min(x, na.rm = T)  )/(max(x, na.rm = T) - min(x, na.rm = T))
}


relative <- function(x){
  if(class(x) != 'numeric'){
    message('non-numeric data provided')
    message('attempting to force numeric')
    x=fn(x)
  }

  x / mean(x,na.rm=T)

}
