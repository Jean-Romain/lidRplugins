#' Check if flightlines are aligned
#'
#' Check if flightlines are well aligned by comparing the average difference between the lowest points
#' of overlapping parts of the flightlines
#'
#' @param las a LAS point cloud
#' @export
check_flightlines_z_alignment <- function(las, res = 3)
{
  PointSourceID <- .N <- NULL

  lay <- lidR:::raster_layout(las, res)
  lay <- lidR:::raster_materialize(lay)
  ids <- unique(las$PointSourceID)
  res <- vector("list", 0)

  for (i in ids) {
    form <- eval(parse(text = paste0("~PointSourceID == ", i)))
    res[[as.character(i)]] <- lidR::pixel_metrics(las, ~min(Z), lay, filter = form)
  }

  n <- las@data[, .N, by = PointSourceID]
  k <- which.max(n$N)
  n$offset <- 0

  for (i in (1:length(ids))[-k]) {
    n$offset[i] = round(mean((res[[i]] - res[[k]])[], na.rm = TRUE), 3)
  }

  n$N <- NULL
  return(n)
}