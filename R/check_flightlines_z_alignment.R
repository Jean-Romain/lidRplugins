#' Check if flightlines are aligned
#'
#' Check if flightlines are well aligned by comparing the average difference between the lowest points
#' of overlapping parts of the flightlines
#'
#' @param las a LAS point cloud
#' @return An adjacency matrix with the offsets between overlapping flightlines
#'
#' @examples
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las = readLAS(LASfile)
#' las <- retrieve_flightlines(las)
#' las$PointSourceID  = las$flightlineID
#' flightlines_z_misalignment_matrix(las)
#' @export
flightlines_z_misalignment_matrix <- function(las, res = 5)
{
  PointSourceID <- NULL

  lay <- lidR:::raster_layout(las, res)
  lay <- lidR:::raster_materialize(lay)
  ids <- unique(las$PointSourceID)
  res <- vector("list", 0)
  data.table::setindex(las@data, PointSourceID)

  if (length(ids) == 1) {
    M <- matrix(0,1,1)
    rownames(M) = ids
    colnames(M) = ids
    return(M)
  }

  for (i in ids) {
    psi <- las@data[PointSourceID == i, .(X,Y,Z)]
    psi <- lidR::LAS(psi, las@header, check = FALSE)
    res[[as.character(i)]] <- lidR:::rasterize_fast(psi, lay, method = "min")
  }

  #plot(terra::rast(res))

  M <- matrix(0, length(ids), length(ids))
  rownames(M) = ids
  colnames(M) = ids

  for (i in 1:length(ids)) {
    for (j in 1:length(ids)) {
      if (i != j) {
        R <- res[[i]] - res[[j]]
        M[i,j] <- round(mean(R[], na.rm = TRUE), 3)
      }
    }
  }

  M[is.nan(M)] <- NA
  return(M)
}


flightlines_z_realignment <- function(las, M)
{
  .N <- PointSourceID <- NULL
  M <- fill_misalignment_matrix(M)
  u <- las@data[, .N, by = PointSourceID]
  idx <- match(las$PointSourceID, u$PointSourceID)
  k <- which.max(u$N)
  offsets <- as.numeric(M[,k])
  newZ <- las$Z - offsets[idx]
  lidR:::quantize(newZ, las[["Z scale factor"]], las[["Z offset"]], by_reference = TRUE)
  las@data[["Z"]] <- newZ
  return(las)
}

fill_misalignment_matrix <- function(M)
{
  Tr <- M
  Tr[is.na(Tr)] <- 0
  Tr[Tr != 0] <- 1
  G <- igraph::graph_from_adjacency_matrix(Tr)

  E <- M
  E[] <- 0

  for (i in 1:nrow(E)) {
    for (j in 1:ncol(E)) {
      if (i != j) {
        paths <- igraph::all_simple_paths(G, i, j, cutoff = 3)
        delta <- vector("list", length(paths))
        E[i,j] <- round(mean(get_all_path_errors(M, paths)), 3)
      }
    }
  }

  return(E)
}

get_all_path_errors = function(M, paths)
{
  deltas <- vector("numeric", length(paths))

  for(i in seq_along(paths)) {
    idx <- as.numeric(paths[[i]])
    e <- 0
    for (j in 1:(length(idx)-1))
      e <- e +  M[idx[j],idx[j+1]]

    deltas[i] <- e
  }

  return(deltas)
}
