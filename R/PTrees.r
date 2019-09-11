#' Individual Tree Detection and Segmentation Algorithm
#'
#' This function is made to be used in \link[lidR:tree_detection]{tree_detection} or \link[lidR:lastrees]{lastrees}.
#' It implements the PTrees algorithm for tree detection and tree segmentation based Vega et al. (2014) (see references).
#' When used in the function \link[lidR:tree_detection]{tree_detection} it runs only the fisrt part of the method i.e. the
#' detection of the trees. When used in  \link[lidR:lastrees]{lastrees} it performs the  whole segmentation (see details).
#'
#' This function has been written by the \code{lidR} authors from the original article. We made our
#' best to implement as far as possible exactly what is written in the original paper but we cannot
#' states that it is the exact original algorithm. Also, minor variations were introduced to fix some
#' issues that were not adressed in the original paper:
#' \itemize{
#' \item Addition of the parameter \code{hmin}: to reduce oversegmentation we introduced a minium height
#' threshold. Points below this thresold cannot initiate new trees during the tree detection and cannot
#' incrase a crown hull during the segmentation.
#' \item Addition of the parameter \code{nmax}: in the original article page 103 figures 5, the number
#' of possible combination is 2^n-n-1. This exponential number of combinations lead, in some cases
#' to an infinite computation time. During the developpement we got cases where the number of combinations
#' to consider was beyond a billion. If the number of tree to consider between two scales is greater
#' than \code{nmax} (i.e. the number of combination is greater than 2^nmax-nmax-1) then the "TreeSegment"
#' from the biggest scale is retained anyway, the smallest scale is considered as dummy.
#' }
#' Notice that to use the PTree strictly as originally described, the point cloud should
#' not be normalized (see reference). In that case the parameter '\code{hmin}' is miningless and can
#' be set to \code{-Inf} for example.
#'
#' @param k integer vector. A serie of k-nearest neighbors to use. In this original paper a k refers
#' to a 'scale' of analyse (see reference).
#'
#' @param hmin scalar. This is an addition from the original paper to limit oversegmentation.
#' Point below this threshold cannot initiate new trees or increase a hull (see details). Set to \code{-Inf}
#' to strictly respect original paper.
#'
#' @param nmax integer. This is an addition from the original paper to protect against uncomputable
#' cases (see details). Set to \code{+Inf} to strictly respect the original paper (not recommended)
#'
#'@author Jasmin Siefert and Jean-Romain Roussel
#'
#' @export
#'
#' @family individual tree segmentation algorithms
#' @family individual tree detection algorithms
#' @family point-cloud based tree segmentation algorithms
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las = readLAS(LASfile, select = "xyz")
#'
#' k = c(30,15)
#' ttops = tree_detection(las, ptrees(k))
#' las   = lastrees(las, ptrees(k))
ptrees = function(k, hmin = 2, nmax = 7L)
{
  assertive::assert_is_numeric(k)
  assertive::assert_all_are_positive(k)
  assertive::assert_all_are_whole_numbers(k)
  assertive::assert_is_a_number(nmax)
  assertive::assert_all_are_whole_numbers(nmax)

  f = function(las)
  {
    context <- tryCatch({get("lidR.context", envir = parent.frame())}, error = function(e) {return(NULL)})
    lidR:::assert_is_valid_context(c("tree_detection", "lastrees"), "ptrees")

    . <- X <- Y <- Z <- treeID <- NULL

    segmentation = context == "lastrees"

    TreeSegments = C_lastrees_ptrees(las, k, hmin, nmax, segmentation)

    if (!segmentation)
    {
      apices = TreeSegments$Apices
      apices = data.table::as.data.table(apices)
      data.table::setnames(apices, names(apices), c("X", "Y", "Z"))
      apices[, treeID := 1:.N]

      output = sp::SpatialPointsDataFrame(apices[, .(X,Y)], apices[, .(treeID, Z)])
      output@proj4string = las@proj4string
      output@bbox = sp::bbox(las)
      return(output)
    }
    else
    {
      return(TreeSegments$treeID)
    }
  }

  class(f) <- c("function", "PointCloudBased", "IndividualTreeSegmentation", "IndividualTreeDetection", "Algorithm", "lidR")
  return(f)
}
