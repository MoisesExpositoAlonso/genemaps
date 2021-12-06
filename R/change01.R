change01<-function (X)
{
    X[X == 1] <- 0.5
    X[X == 2] <- 1
    X
}

change11neg<-function (X){
  X[X == 1] <- 0
  X[X == 0] <- (-1)
  X[X == 2] <- 1
  X
}
