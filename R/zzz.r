#' @useDynLib lidRplugins
#' @importFrom Rcpp sourceCpp
#' @import data.table
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("lidRplugins", libpath)
}