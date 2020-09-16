#' Individual Tree Detection Algorithm
#'
#' This function is made to be used in \link[lidR:find_trees]{find_trees}. It implements an algorithms for tree
#' detection based on a method described in Eysn et al (2015) (see references) and propably proposed
#' originaly by someone else (we did not find original publication). This is a local maximum filter
#' applied on a multi-canopy height model (see details).\cr\cr
#' Notice: tree tops returned are the true highest points within a given pixel whenever the CHMs where
#' computed with the 95th percentile of height. Otherwise these maxima are not true maxima and cannot
#' be used in subsequent segmentation algorithms.
#'
#' Description adapted from Eysn et al (2015), page 1728, section 3.1.3 Method #3\cr\cr
#' The method is based on iterative canopy height model generation (CHM) and local maximum filter (LMF)
#' detection within a moving window for various CHMs. The method works in two general steps, which are
#' (a) sequential identification of potential trees and (b) filtering of the extracted potential trees.\cr
#' \itemize{
#' \item Step (a): From the normalized point cloud, an initial CHM is created by assigning the 95th
#' height percentile within each raster cell. Based on this CHM, LM are detected and the found
#' positions and heights are stored in a database. For the next iteration, points in the uppermost layer
#' of the normalized ALS data are eliminated. The “eliminating” layer is defined as a band below the
#' current CHM. Based on the filtered data, a new CHM is created, LM are extracted, and the
#' LM parameters are added to the database. This procedure is carried out sequentially until all points
#' are removed from the normalized point cloud.
#' \item Step (b): All detected LM in the database are sorted by decreasing heights.
#' The highest LM is considered a detected tree. For each following LM, the LM is considered a
#' detected tree if there is no detected tree within a given 2D distance as well as a given 3D distance.
#' }

#' @param res numeric. Resolution of the CHM based on the 95th percentile
#'
#' @param layer_thickness numeric. The “eliminating” layer is defined as a band of \code{layer_thickness} m
#' below the current CHM (see details).
#'
#' @param dist_2d numeric. 2D distance threshold. A local maximum is considered a detected tree
#' if there is no detected tree within this 2D distance (see details).
#'
#' @param dist_3d numeric. 3D distance threshold. A local maximum is considered a detected tree
#' if there is no detected tree within this 3D distance (see details).
#'
#' @param use_max logical. The CHMs are computed with the 95th percentiles of height in the original
#' description. If \code{use_max = TRUE} it uses the 100th percentiles (max height) and thus does
#' not implies any sorting algorithm. The algoithm is therefore 5 to 10 times faster.
#'
#' @param ... supplementary parameters to be passed to \link[lidR:lmf]{lmf} that is used internally
#' to find the local maxima.
#'
#' @export
#'
#' @references
#' Eysn, L., Hollaus, M., Lindberg, E., Berger, F., Monnet, J. M., Dalponte, M., … Pfeifer, N. (2015).
#' A benchmark of lidar-based single tree detection methods using heterogeneous forest data from the
#' Alpine Space. Forests, 6(5), 1721–1747. https://doi.org/10.3390/f6051721
#'
#' @family individual tree detection algorithms
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las = readLAS(LASfile)
#'
#' ttops = find_trees(las, multichm(res = 1, ws = 5))
#'
#' x = plot(las)
#' add_treetops3d(x, ttops)
multichm = function(res = 1, layer_thickness = 0.5, dist_2d = 3, dist_3d = 5, use_max = FALSE, ...)
{
  lidR:::assert_is_a_number(res)
  lidR:::assert_is_a_number(layer_thickness)
  lidR:::assert_is_a_number(dist_2d)
  lidR:::assert_is_a_number(dist_3d)
  lidR:::assert_all_are_positive(res)
  lidR:::assert_all_are_positive(layer_thickness)
  lidR:::assert_all_are_positive(dist_2d)
  lidR:::assert_all_are_positive(dist_3d)

  f = function(las)
  {
    context <- tryCatch({get("lidR.context", envir = parent.frame())}, error = function(e) {return(NULL)})
    lidR:::assert_is_valid_context(lidR:::LIDRCONTEXTITD, "multichm")

    . <- X <- Y <- Z <- treeID <- NULL

    dist_2d <- dist_2d^2
    dist_3d <- dist_3d^2

    las_copy <- lidR::LAS(las@data[, .(X,Y,Z)], las@header)
    LM       <- list()
    chm      <- lidR::grid_canopy(las, res, lidR::p2r())
    i        <- 1
    p        <- list(...)
    hmin     <- if (is.null(p$hmin)) formals(lidR::lmf)$hmin else p$hmin

    while (!lidR::is.empty(las_copy))
    {
      if (use_max)
        chm95 <- lidR::grid_canopy(las_copy, res, lidR::p2r())
      else
        chm95 <- lidR::grid_metrics(las_copy, ~stats::quantile(Z, probs = 0.95), res)

      if (max(chm95[], na.rm = TRUE) > hmin)
      {
        lm       <- lidR::find_trees(chm95, lidR::lmf(...))
        colnames(lm@coords) <- c("X", "Y")
        lm       <- raster::as.data.frame(lm)
        data.table::setDT(lm)
        LM[[i]]  <- lm
        las_copy <- lidR::lasmergespatial(las_copy, chm95, "chm95")
        las_copy <- lidR::lasfilter(las_copy, Z < chm95 - layer_thickness)

        i <- i + 1
      }
      else
        las_copy <- methods::new("LAS")
    }

    LM <- data.table::rbindlist(LM)
    data.table::setorder(LM, -Z)
    LM <- unique(LM, by = c("X", "Y"))

    detected = logical(nrow(LM))
    detected[1] = TRUE

    for (i in 2:nrow(LM))
    {
      distance2D = (LM$X[i] - LM$X[detected])^2 + (LM$Y[i] - LM$Y[detected])^2
      distance3D = distance2D + (LM$Z[i] - LM$Z[detected])^2

      if (!any(distance2D < dist_2d) & !any(distance3D < dist_3d))
      {
        detected[i] = TRUE
      }
    }

    detected = LM[detected]
    detected[, treeID := 1:.N]

    output <- sp::SpatialPointsDataFrame(detected[, 3:4], detected[, 1:2], proj4string = las@proj4string)
    return(output)
  }

  class(f) <- c(lidR:::LIDRALGORITHMITD, lidR:::LIDRALGORITHMPOINTCLOUDBASED)
  return(f)
}
