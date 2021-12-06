meancent.mat<-function (X.) 
{
    X = apply(X., 2, function(x) {
        mu = mean(x)
        x - mu
    })
    X
}
meanvarcent<-function (x) 
{
    if (class(x) != "numeric") {
        message("non-numeric data provided")
        message("attempting to force numeric")
        x = fn(x)
    }
    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
}

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