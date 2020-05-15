#' Individual Tree Detection Algorithm
#'
#' This function is made to be used in \link[lidR:find_trees]{find_trees}. It implements an
#' experimental algorithms for tree detection based on a several ideas from the litterature. First it
#' select the highest points in each cell of a 1 m grid to reduce the amount of data and considerably
#' improve speed, then it performs a local maximum filter to find tree tops. To finish it applies the
#' filtering rule from \link{multichm} (step (b))
#'
#' @param ws numeric or function. Length or diameter of the moving window used to the detect the local
#' maxima in the unit of the input data (usually meters). If it is numeric a fixed windows size is used.
#' If it is a function, the function determines the size of the window at any given location on the canopy.
#' The function should take the height of a given pixel or points as its only argument and return the
#' desired size of the search window when centered on that pixel/point.
#'
#' @param hmin numeric. Minimum height of a tree. Threshold below which a pixel or a point
#' cannot be a local maxima. Default 2.
#'
#' @param dist_2d numeric. 2D distance threshold. A local maximum is considered a detected tree
#' if there is no detected tree within this 2D distance.
#'
#' @export
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las = readLAS(LASfile)
#'
#' ttops = tree_detection(las, lmfx(ws = 3))
#'
#' x = plot(las)
#' add_treetops3d(x, ttops)
lmfx = function(ws, hmin = 2, dist_2d = 3)
{
  f = function(las)
  {
    context <- tryCatch({get("lidR.context", envir = parent.frame())}, error = function(e) {return(NULL)})
    lidR:::assert_is_valid_context(lidR:::LIDRCONTEXTITD, name = "lmfx")

    . <- X <- Y <- Z <- treeID <- NULL

    dist_2d = dist_2d^2

    if (assertive::is_a_number(ws))
    {
      # nothing to do
    }
    else if (assertive::is_function(ws))
    {
      n  = nrow(las@data)
      ws = ws(las@data$Z)

      if (!is.numeric(ws)) stop("The function 'ws' did not return correct output. ", call. = FALSE)
      if (any(ws <= 0))    stop("The function 'ws' returned negative or nul values.", call. = FALSE)
      if (anyNA(ws))       stop("The function 'ws' returned NA values.", call. = FALSE)
      if (length(ws) != n) stop("The function 'ws' did not return correct output.", call. = FALSE)
    }
    else
      stop("'ws' must be a number or a function", call. = FALSE)

    . <- X <- Y <- Z <- treeID <- NULL
    las = lidR::decimate_points(las, lidR::highest(1))
    is_maxima = lidR:::C_lmf(las, ws, hmin, TRUE, lidR:::getThread())
    LM = las@data[is_maxima, .(X,Y,Z)]

    data.table::setorder(LM, -Z)

    detected = logical(nrow(LM))
    detected[1] = TRUE

    for (i in 2:nrow(LM))
    {
      distance2D = (LM$X[i] - LM$X[detected])^2 + (LM$Y[i] - LM$Y[detected])^2
      #distance3D = distance2D + (LM$Z[i] - LM$Z[detected])^2

      if (!any(distance2D < dist_2d))# & !any(distance3D < dist_3d))
      {
        detected[i] = TRUE
      }
    }

    detected = LM[detected]
    detected[, treeID := 1:.N]

    output = sp::SpatialPointsDataFrame(detected[, 1:2], detected[, 3:4])
    output@proj4string = las@proj4string
    output@bbox = sp::bbox(las)
    return(output)
  }

  class(f) <- c(lidR:::LIDRALGORITHMITD)
  return(f)
}
