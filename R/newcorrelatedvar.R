

newcorrelatedvar <- function(Y,rho,samevariance=TRUE) {
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
  if(samevariance) x=x*sd(Y)
  return(x)
}

addnoise<-function(Y,noise=0.5,bias=0){
  Y + rnorm(length(Y),mean(Y)*bias,sd = sd(Y))*noise
}
