#' Test if points are part of a neighborhood that is approximately coplanar
#'
#' An approximate coplanarity test. For each point it looks for the k-nearest neighbors. It computes the
#' eigenvalues of the covariance matrix. The points that meet the following criteria are labeled as
#' approximately coplanar:\cr
#' \deqn{a_2 > (th_1*a_1) && (th_2*a_2) > a_3}
#' \eqn{a_1, a_2, a_3} are the eigenvalues of a neighborhood of points (defined by k-nearest neighbors)
#' in ascending order.\cr\cr
#' If \code{th2 = 0} the function turns to be a colinearity test. The points that meet the following criteria are labeled as
#' approximately colinear:\cr
#' \deqn{th_1*a_2 < a_3 && (th_1*a_1) < a_3}
#'
#' @return A LAS object with a new column names \code{Coplanar} (or \code{Colinear} if \code{th2 = 0})
#' that indicates those points that are part of a neighborhood that is approximately coplanar/colinear
#'  (TRUE) or not (FALSE).
#'
#' @param las an object of class LAS
#' @param k interger. The number of k-nearest neighbors.
#' @param th1 numeric. The threshold to be applied to the smallest eigenvalue.
#' @param th2 numeric. The threshold to be applied to the second smallest eigenvalue.
#'
#' @export
lascoplanar = function(las, k = 8, th1 = 25, th2 = 6)
{
  stopifnot(nrow(las@data) > k, th1 > 0, th2 >= 0)
  name <- if (th2 == 0) "Colinear" else "Coplanar"
  out <- C_lascoplanar(las, k, th1, th2)
  las@data[[name]] <- out
  return(las)
}