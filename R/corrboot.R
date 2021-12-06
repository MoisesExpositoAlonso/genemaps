

corrboot<-function(dat=mtcars[, c(1, 3)], R=100,method='s'){
dat<-data.frame(as.matrix(na.omit(dat)))
N <- nrow(dat)

cor.orig <- cor(dat)[1,2]
cor.boot <- NULL

for (i in 1:R) {
  idx <- sample.int(N, N, replace = TRUE) 
  cor.boot[i] <- cor(dat[idx, ],method=method)[1,2] 
}
return(list(mean=mean(cor.boot),
            lower=quantile(cor.boot,probs = 0.025) ,
            upper=quantile(cor.boot,probs = 0.975)
            ))
}
