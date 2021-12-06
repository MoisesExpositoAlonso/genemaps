#' #' Run individual-based simulations
#' #'
#' #' @param types
#' #' @param proportions
#' #' @param fitness
#' #' @param t
#' #' @param d
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' populationsimC <-  function(
#'                             types,
#'                             individuals,
#'                             fitness,
#'                             t,
#'                             d,
#'                             verbose='false') {
#'
#'
#'     result=.Call('gws_populationsimC', PACKAGE = 'gws',
#'                   types=as.numeric(types),
#'                   individuals=as.numeric(individuals),
#'                   fitness=as.numeric(fitness),
#'                   t=as.numeric(t),
#'                   d=as.numeric(d),
#'                   verbose=as.logical(verbose)
#'                   )
#'     rownames(result)<-types
#'     return(result)
#' }
#'
#'
#' #' Run individual-based simulations in a grid search for the drift parameter
#' #'
#' #' @param types
#' #' @param individuals
#' #' @param fitness
#' #' @param t
#' #' @param dmin
#' #' @param dmax
#' #' @param dstep
#' #' @param verbose
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' populationsimCgrid <- function(types, individuals, fitness, t = 20, dmin = 0, dmax = 0.1, dstep = 0.0001, verbose = FALSE) {
#'     .Call('gws_populationsimCgrid', PACKAGE = 'gws', types, individuals, fitness, t, dmin, dmax, dstep, verbose)
#'
#'       result=.Call('gws_populationsimCgrid', PACKAGE = 'gws',
#'                   types=as.numeric(types),
#'                   individuals=as.numeric(individuals),
#'                   fitness=as.numeric(fitness),
#'                   t=as.numeric(t),
#'                   dmin=as.numeric(dmin),
#'                   dmax=as.numeric(dmax),
#'                   dstep=as.numeric(dstep),
#'                   verbose=as.logical(verbose)
#'                   )
#'     # rownames(result)<-types
#'     return(result)
#'
#' }
