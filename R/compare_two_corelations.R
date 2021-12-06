# fishertransformation<-function(r){
#   # Transformation of r to z by Fisher
#   z= .5*(log(1+r)-log(1-r))
#   return(z)
# }
# fishertransformation(0.5)
# 
# N1=50
# N2=50
# r2=0.6
# r1=0.4
# test_two_correlations<-function(r1,r2=0,N1,N2=N1){
#   # Test the difference between two r coefficients based on Z scores
#   z1=fishertransformation(r1)
#   z2=fishertransformation(r2)
#   Zobserved = (z1-z2) / ( sqrt( (1 / N1-3) + (1/N2-3) ) )
#   pnorm(Zobserved)
# }
