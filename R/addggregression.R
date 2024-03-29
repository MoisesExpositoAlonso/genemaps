#' Add a regression or loess to an existing ggplot
#'
#' @param myplot The ggplot object that wants to be annotated with a regression
#'
#' @return Returns the pots with the additions
#'
#' @seealso ggdotscolor
#'
#' @examples
#' x=rnomr(1,1,100)
#' y=rnorm(1,2,100)
#' x=rnorm(1,3,100)
#' ggdotscolor(x=x,y=y,colorvar=z) %>% addggregression()
#'
#' @export
#'


addggregression<-function(myplot  ,doregression=T ,doloess=F, docorrelation=F, span=if(doloess==T){ c(.5,1,1.5,2)} ,colour="black", ...) {

  if(doregression==T & docorrelation==F){
    myplot<- myplot+ stat_smooth (method="glm",colour=colour, ...) +
      # annotate("text",  x=Inf, y = Inf, label = lm_eq(myplot$data$y,myplot$data$x), vjust=1, hjust=1,parse=TRUE)
      ggtitle(TeX(lm_eq(myplot$data$y,myplot$data$x,tex = T)) )
  }
  if(doloess==T & docorrelation==F){
  for(span in span){
      myplot<-  myplot + stat_smooth(method="loess",span=span, colour=colour, ...)
    }
  }
	if(docorrelation==T){
	    myplot<- myplot+ stat_smooth (method="glm",colour=colour, ...) +
	    	ggtitle(TeX(r2_eq(myplot$data$y,myplot$data$x,tex = T)) )

	}
	
  return(myplot)
}
