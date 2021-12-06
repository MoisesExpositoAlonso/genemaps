
qqGWA<-function(pvals,dir=getwd(),name="",save=F){

  observed <- sort(pvals)
  lobs <- -(log10(observed))

  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))

  if(save==T){
    pdf(paste("qqplot",name,".pdf",sep=""), width=6, height=6)
  }

  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,max(lobs)+1), ylim=c(0,max(lobs)+1), las=1, xaxs="i", yaxs="i", bty="l")
  points(lexp, lobs, pch=23, cex=.4, bg="black")
  abline(lm(lobs~lexp ),col="black",lty="dashed")
  title( paste("lambda",round(coefficients(lm(lobs~lexp ))[2],digits=3) ) )
  if(save==T){
    message(paste("qqplot printed to: ",dir))
    dev.off()
}}
