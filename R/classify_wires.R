#' Classify wire
#'
#' Classify wires using the output of \link{track_wires}. Using wire profiles and some knowledge
#' on the geometry of the towers it creates a cloth below the wires and classifies everything that
#' is above this cloth will be classifies a wire exept the transmission towers. This the transmission
#' tower classification must be performed first.
#'
#' @param las An object of class LAS.
#' @param wires A \code{SpatialPointsDataFrame} returned by \link{track_wires}
#' @param dtm A RasterLayer. The digital terrain model is need to get the relative elevations
#'
#' @examples
#' \donttest{
#' # A simple file with wires already clipped from 4 files + shapefile
#' # of the network
#' LASfile <- system.file("extdata", "wires.laz", package="lidRplugins")
#' wireshp <- system.file("extdata", "wires.shp", package="lidRplugins")
#' las <- readLAS(LASfile, select = "xyzc")
#' network <- shapefile(wireshp)
#'
#' # Ugly dtm because data is an Y
#' dtm <- grid_terrain(las, 2, tin())
#'
#' towers <- find_transmissiontowers(las, network, dtm, "waist-type")
#' las <- classify_transmissiontowers(las, towers, dtm)
#' wires <- track_wires(towers, network, dtm, "waist-type")
#' las <- classify_wires(las, wires, dtm)
#'
#' plot(las, color = "Classification")
#' }
#' @family electrical network
#' @export
classify_wires = function(las, wires, dtm)
{
  UseMethod("classify_wires", las)
}

#' @export
classify_wires.LAS = function(las, wires, dtm)
{
  classify_from_virtual = FALSE

  SECTIONS = unique(wires$section)

  las2 <- lasmergespatial(las, dtm, "dtm")
  las2@data$ID <- 1:npoints(las)

  for (section in SECTIONS)
  {
    wire = wires[wires$section == section,]
    tower.spec = get_tower_spec(wire$type[1])
    wire$type <- NULL

    thresholds = 0
    if (tower.spec$wire.layers == 1) {
      thresholds = 6
    } else {
      thresholds = (tower.spec$wire.layers - 1) * tower.spec$wire.distance + 5
    }

    lwires <- sp::SpatialLines(list(sp::Lines(list(sp::Line(wire@coords)), ID = "1")))
    pwires <- rgeos::gBuffer(lwires, width = 0.5*tower.spec$length[2], capStyle = "SQUARE")

    sub <- lasclip(las2, extent(pwires))
    layout <- lidR:::rOverlay(sub, 10)
    cloth <- raster::rasterize(wire, layout)$z

    ker <- matrix(1,3,3)
    for (k in 1:2)
      cloth <- raster::focal(cloth, ker, fun = stats::median, na.rm = TRUE, pad = T)

    sub <- lasmergespatial(sub, pwires, "pwires")
    sub <- lasmergespatial(sub, cloth, "cloth")
    sub$Classification[sub$id == 1 & sub$Z > sub$cloth - thresholds & sub$Classification != lidR::LASTRANSMISSIONTOWER] <- lidR::LASWIRECONDUCTOR
    ids = sub$ID[sub$Classification == lidR::LASWIRECONDUCTOR]
    las@data[["Classification"]][ids] <- lidR::LASWIRECONDUCTOR

    # Check if some wire are classified from virtual tracks
    if (any(wire$virtual == 1))
    {
      # Some points are classified and no buffer in the las
      if (length(ids) > 0 && is.null(las@data[["buffer"]]))
      {
        classify_from_virtual = TRUE
      }
      # Some points are classified and a buffer in the las
      else if (length(ids) > 0 && !is.null(las@data[["buffer"]]))
      {
        # Somme classified points are not from the buffer
        if (!all(las@data[["buffer"]][ids] == 0))
          classify_from_virtual = TRUE
      }
    }
  }

  if (classify_from_virtual)
    warning("Some points were classified as wires using virtual tracks and are likely to be poorly classified", call. = FALSE)

  return(las)
}