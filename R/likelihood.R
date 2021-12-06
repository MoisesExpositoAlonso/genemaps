RGamma_lik<-function(a,b,y,minc){
  L = (a-1) * log(y) - y*(1/b) - log(gamma(a+minc)) -  a * log(1/b);
  return(L)
}


ld.Cauchy.lik<-function (theta, y)
{
    e = theta[1]
    sel = theta[-1]
    ldexpected <- ldCexpect(sel, e) %>% cleanld()
    location <- median(ldexpected)
    scale <- IQR(ldexpected)/2
    L <- suppressWarnings(sum(dcauchy(x, location, scale, log = T)))
    return(L)
}


ld.exp.lik<-function (theta, y)
{
    e = theta[1]
    sel = theta[-1]
    ldexpected <- lowermat2vec(ldCexpect(s, e))
    rate = 1/mean(ldexpected)
    return(sum(-dexp(x = x, rate = theta, log = T)))
}

ld.Gauss.lik<-function (theta, y)
{
    e = theta[1]
    sel = theta[-1]
    ldexpected <- lowermat2vec(ldCexpect(sel, e))
    mu <- mean(ldexpected)
    sigma2 <- var(ldexpected)
    n <- length(y)
    logl <- -0.5 * n * log(2 * pi) - 0.5 * n * log(sigma2) -
        (1/(2 * sigma2)) * sum((y - mu)^2)
    return(-logl)
}
