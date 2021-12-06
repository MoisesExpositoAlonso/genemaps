#' Title
#'
#' @param mod model of type MCMCglmm. It needs to be runed with random effects and the flag pl=TRUE, pr=TRUE
#' @param listrandom strings of the name used for the random factor of interest
#'
#' @return
#' @export
#'
#' @examples

extract_MCMCglmm_randomeffects<-function(mod, listrandom){

listrandom=paste0(listrandom,'.')

alleffects=colnames(mod$Sol)

solutions<- data.frame(effect=alleffects, posteriormode=apply(mod$Sol, 2, posterior.mode))

targeteffects<-alleffects[
  unlist(
  lapply(listrandom, function(x) grep(x,alleffects))
)]

mysol<-filter(solutions, effect %in% targeteffects )

return(mysol)
}
