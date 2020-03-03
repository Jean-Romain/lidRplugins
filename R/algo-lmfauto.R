#' Individual Tree Detection Algorithm
#'
#' This function is made to be used in \link{find_trees}. It implements a fast and parameter-free
#' algorithm for individual tree detection for broad coverage. It is based on two local maximum filters
#' (LMF). The first pass performs a very rough estimation of the number of trees with a fixed window
#' size. Based on this rough estimate it automatically computes a variable windows size LMF with workable
#' parameters. This algorithm is made to process wide areas rather than small plots. See references
#' for more details.
#'
#' @param plot logical set it to \code{TRUE} if processing a plot instead of a large area. What changes
#' is the estimation of the local number of trees. It should be based on the local neighborhood for the general
#' case but this does not make sense for a plot.
#' @param hmin numeric. Minimum height of a tree. Threshold below which a point cannot be a local
#' maxima. Default is 2.
#'
#' @references Jean-Romain Roussel, Francesco Pirotti, Luiz Carlos Estraviz Rodriguez,
#' Jean-François Bourdon, Antoine Lebœuf, Marc-Olivier Lemonde, Alexis Achim. Development of an
#' auto-adaptive individual tree detection algorithm for airborne LiDAR data (in prep.)
#'
#' @family individual tree detection algorithms
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#' ttops <- find_trees(las, lmfauto())
#'
#' x = plot(las)
#' add_treetops3d(x, ttops)
#' @export
lmfauto = function(plot = FALSE, hmin = 2)
{
  f = function(las)
  {
    lidR:::assert_is_valid_context(lidR:::LIDRCONTEXTITD, "lmfauto")

    # Step 1: detection with a fixed 5 m windows size
    ttop5 <- lidR::find_trees(las, lidR::lmf(5))

    # Step 2: raw/rough/poor estimate of number of trees per ha in
    # the local neighourhood

    if (plot)
    {
      # Limit case if we are not processing a wide area
      A     <- lidR::area(las)
      d     <- nrow(las@data)/A
      Aha   <- 10000/A
      ntop5 <- nrow(ttop5)*Aha
    }
    else
    {
      # The real algorithm
      A     <- 400
      Aha   <- 10000/A
      x     <- ttop5@coords[,1]
      y     <- ttop5@coords[,2]
      ntop5 <- lidR:::C_count_in_disc(x, y, las@data$X, las@data$Y, sqrt(A/pi), lidR:::getThread())
      ntop5 <- ntop5*Aha
    }

    # Step 3: estimate the window size of a variable window size LMF as a function
    # of the number of trees in the local neighborhood.
    . <- X <- Y <- Z <- treeID <- NULL

    ws <- lmfauto_ws(las@data$Z, ntop5)
    lm <- lidR:::C_lmf(las, ws, hmin, TRUE, lidR:::getThread())
    return(lm)
  }

  class(f) <- c("PointCloudBased", "IndividualTreeDetection", "Algorithm", "lidR")
  return(f)
}

lmfauto_ws = function(x, n, d = 10)
{
  s <- length(n)
  above200 <- n > 200
  above300 <- n > 300

  a <- rep(3.5, s)
  b <- rep(4, s)
  a[above200] <- 2.5
  b[above200] <- 3.5
  a[above300] <- 1.5
  b[above300] <- 2.5

  if (d < 4)
  {
    a <- a + 1.25
    b <- b + 1.25
    a[above200] <- a[above200] - 0.5
    b[above200] <- b[above200] - 0.5
    a[above300] <- a[above300] - 0.25
    b[above300] <- b[above300] - 0.25
  }

  llim  <- 2
  ulim  <- 20
  slope <- (b - a)/(ulim - llim)
  intercept <- a - 2*slope
  ws <- slope*x + intercept
  ws[x < llim] <- a[x < llim]
  ws[x < llim] <- b[x < llim]
  return(ws)
}
