#' Classify the transmission towers
#'
#' Attribute the class 15 to points that correspond to transmission towers using the positionning
#' of the towers given by \link{find_transmissiontowers}. Using the transmission towers types and their
#' orientation it computes the spatial box that emcompasses the towers and classify points within this
#' rectangle as 'transmission tower'.
#'
#' @param las An object of class LAS
#' @param towers SpatialPointsDataFrame returned by \link{find_transmissiontowers}.
#' @param dtm A RasterLayer. The digital terrain model is useful to find the bottom of the towers
#' @param threshold numeric. Height above ground. Points below this elevation are not classified
#' as transmission towers.
#'
#' @return A LAS object with an updated classification.
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
#'
#' plot(las, color = "Classification")
#' }
#' @family electrical network
#' @export
classify_transmissiontowers = function(las, towers, dtm, threshold = 2)
{
  UseMethod("classify_transmissiontowers", las)
}

#' @export
classify_transmissiontowers.LAS = function(las, towers, dtm, threshold = 2)
{
  towers <- tower.boundingbox(towers)
  tmp <- lidR::lasmergespatial(las, towers, "towers")
  tmp <- lidR::lasnormalize(tmp, dtm)
  las@data$Classification[tmp$Z > threshold & tmp$towers == TRUE] <- lidR::LASTRANSMISSIONTOWER
  return(las)
}

tower.boundingbox = function(towers)
{
  dtm <- NULL

  if (length(towers) == 0L)
  {
    data = data.frame(maxZ = numeric(0), minZ = numeric(0), deflection = integer(0))
    out = sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list()), data)
    raster::projection(out) <- raster::projection(towers)
    out@bbox = towers@bbox
    return(out)
  }

  lines <- vector("list", length(towers))
  for (i in 1:length(towers))
  {
    tower <- towers[i,]
    tower.spec <- get_tower_spec(tower$type)

    height <- tower.spec$width[2]
    width <- tower.spec$length[2]
    hwidth <- width/2
    hheight <- height/2
    p1 <- tower@coords
    ux <- tower$ux
    uy <- tower$uy
    orientation <- matrix(c(ux, uy, uy, -ux), ncol = 2)

    p2 <- p1
    p2[,1] <- p2[,1] + orientation[1,2] * (hwidth - hheight)
    p2[,2] <- p2[,2] + orientation[2,2] * (hwidth - hheight)
    p1[,1] <- p1[,1] - orientation[1,2] * (hwidth - hheight)
    p1[,2] <- p1[,2] - orientation[2,2] * (hwidth - hheight)

    lines[[i]] <- sp::Lines(sp::Line(rbind(p1, p2)), as.character(i))
  }

  tower.orientation <- sp::SpatialLines(lines, proj4string = towers@proj4string)
  #plot(tower.orientation, add = T)

  # Compute an extent for the tower by buffering the lines
  towers.extent <- rgeos::gBuffer(tower.orientation, width = tower.spec$width[2]/2, capStyle = "SQUARE", byid = T)
  towers.extent$maxZ <- towers$Z
  towers.extent$minZ <- sapply(raster::extract(dtm, towers.extent), mean)

  return(towers.extent)
}
