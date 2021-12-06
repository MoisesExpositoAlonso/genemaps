
CIcor<-function(r=0.5,num=50){
  stderr = 1.0 / sqrt(num - 3)
  delta = 1.96 * stderr
  lower = tan(atan(r) - delta)
  upper = tan(atan(r) + delta)
  # print "lower %.6f upper %.6f" % (lower, upper)
  return(cbind(lower,upper))
}
