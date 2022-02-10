#' Track the wires of the powerlines using the tower positions
#'
#' Assuming the coordinates, the elevation and the type of trasmission tower are known the functions
#' tracks the wires of the powerlines. To achieve this task it computes the catenary equation of the
#' wires between two consecutive towers. Striclty speaking the function is not able to find the wires
#' instead it guesses their equation that can be computed deterministically from the tower coordinates,
#' their height and their type.
#'
#' @param towers \code{SpatialPointsDataFrame} containing the positions of the towers as returned by
#' \link{find_transmissiontowers}
#' @param powerline A \code{SpatialLines*} that map the electrical network accurately. The method may
#' be improved later to get rid of this information.
#' @param dtm A RasterLayer. Digital Terrain Model is useful to find a relevant elevation for virtual
#' towers
#' @param type One of "waist-type" or "double-circuit" according to
#' \href{http://www.hydroquebec.com/learning/transport/types-pylones.html}{Hydro-Quebec}. Can also
#' be a list with custom specifications. See examples.
#' @param debug logical. Plot the different steps of the algorithm so one can try to figure out what
#' is going wrong.
#'
#' @return  A \code{SpatialPointsDataFrame} that actually represents 3D lines with several attributes
#' per points. \code{Z} the elevation, \code{virtual} tells if the wire has been found with two
#' consecutive tower and is thus accurate or if only one tower was used and in this case the wire is
#' a pure guess, \code{section} attributes an ID to each wire section i.e. between two towers, \code{ID}
#' attributes and ID to each powerline \code{type} store the transmission tower type.
#'
#' @references
#' Roussel J, Achim A, Auty D. 2021. Classification of high-voltage power line structures in low density
#' ALS data acquired over broad non-urban areas. PeerJ Computer Science 7:e672 https://doi.org/10.7717/peerj-cs.672
#'
#' @examples
#' \donttest{
#' # A simple file with wires already clipped from 4 files + shapefile
#' # of the network
#' LASfile <- system.file("extdata", "wires.laz", package="lidRplugins")
#' wireshp <- system.file("extdata", "wires.shp", package="lidRplugins")
#' dtmtif  <- system.file("extdata", "wire-dtm.tif", package="lidRplugins")
#' las <- readLAS(LASfile, select = "xyzc")
#' network <- sf::st_read(wireshp)
#' dtm <- raster::raster(dtmtif)
#'
#' towers <- find_transmissiontowers(las, network, dtm, "waist-type")
#' wires <- track_wires(towers, network, dtm, "waist-type")
#'
#' col <- c("red", "blue", "forestgreen", "darkorchid", "darkorange")[wires$ID]
#' col[wires$virtual & col == "red"] <- "pink"
#' col[wires$virtual & col == "blue"] <- "lightblue"
#' col[wires$virtual & col == "forestgreen"] <- "lightgreen"
#' col[wires$virtual & col == "darkorchid"] <- "plum"
#' col[wires$virtual & col == "darkorange"] <- "goldenrod1"
#'
#' plot(header(las))
#' plot(towers, add = T, col = towers$deflection + 1)
#' plot(wires, col = col, add = T, cex = 0.1)
#'
#' plot(las, clear_artifacts = FALSE)
#' rgl::points3d(wires@coords[,1], wires@coords[,2], wires$z, col = col, size = 5)
#' }
#' @family electrical network
#' @export
track_wires <- function(towers, powerline, dtm, type = c("waist-type", "double-circuit"), debug = FALSE)
{
  if (is(powerline, "sf") | is(powerline, "sfc")) powerline <- sf::as_Spatial(powerline)

  tower.spec <- get_tower_spec(type)

  tlocation <- towers
  proj <- tlocation@proj4string

  # If 0 tower the question is closed: return nothing
  if (length(tlocation) == 0)
  {
    coords <- matrix(0, ncol = 2)
    data   <- data.frame(z = numeric(1), virtual = integer(1), ID = 0)
    output <- sp::SpatialPointsDataFrame(coords, data, proj4string = tlocation@proj4string)
    return(output[0,])
  }

  if (debug)
  {
    opar = graphics::par("mfrow")
    graphics::par(mfrow = c(2,3))
    on.exit(graphics::par(mfrow = opar))
  }

  # Crop the lines to the extent of the ROI
  pwll <- raster::crop(powerline, extent(dtm))

  if (debug)
  {
    plot(raster::extent(dtm), main = paste0("Raw powerline network"), asp = 1)
    plot(pwll, add = T, col = 1:length(pwll))
  }

  # This starts like the tower detection by fixing the shapefile
  pwll <- gJoinLines(pwll, 2)
  pwll <- rgeos::gSimplify(pwll, 40)
  spwll <- gSplitLines(pwll)
  spwlp <- rgeos::gBuffer(spwll, width = 125, byid = T, capStyle = 'SQUARE')

  if (debug)
  {
    plot(raster::extent(dtm), main = "Post-processed lines", asp = 1)
    #plot(as(spwll, "SpatialPoints"), add = T)
    plot(spwll, add = T, col =  1:length(spwll))
    plot(spwlp, add = T, border =  1:length(spwlp), lty = 3)
  }

  # Decompose the wires with different orientations into linear sections
  if (any(towers$deflection))
  {
    angles <- unique(round(tlocation$theta, 2))
    tlocations <- vector("list", length(angles))
    for (i in 1:length(angles)) tlocations[[i]] <- tlocation[round(tlocation$theta,2) == angles[i],]
  } else {
    tlocations <- list(tlocation)
  }

  if (length(tlocations) != length(spwll))
    stop("Internal error: different number of sections.", call. = FALSE)

  if (debug)
  {
    plot(raster::extent(dtm), main = "Detection of the powerlines", asp = 1)
  }

  # Loop on each section
  ID <- 1                    # ID for each line
  SECTION <- 0               # ID for each section
  HXY <- vector("list", length(tlocations))
  for (kk in 1:length(tlocations))
  {
    # Get the section
    tlocation <- tlocations[[kk]]
    pwlp <- spwlp[kk,]
    #plot(tiles)
    #plot(tlocation, add = T)
    #plot(pwlp, add = T)

    # Initialize vars
    n <- nrow(tlocation)
    posx <- tlocation@coords[,1]
    posy <- tlocation@coords[,2]
    wire <- vector("list", n)
    #plot(tiles)
    #plot(tlocation, add = T, col = tlocation$deflection +1)
    #arrows(posx, posy, posx + 100*tlocation$ux, posy + 100*tlocation$uy, length = 0.05)

    # Draw 1000 m lines passing throught each tower
    # Then buffer to create a polygons that encommpass each wire line
    for (i in 1:n)
    {
      l <- 500
      x <- posx[i]
      y <- posy[i]
      ux <- tlocation$ux[i]
      uy <- tlocation$uy[i]
      coords <- matrix(c(x - l*ux, y - l*uy, x + l*ux, y + l*uy), ncol = 2, byrow = T)
      wire[[i]] <- sp::Lines(list(sp::Line(coords)), ID = as.character(i))
    }

    lwires <- sp::SpatialLines(wire, proj4string = proj)
    crlwires <- raster::crop(lwires, extent(dtm) - 1)
    pwires <- rgeos::gBuffer(lwires, width = 0.3*tower.spec$length[2], capStyle = "SQUARE")
    pwires <- raster::crop(pwires, pwlp)
    pwires <- sp::disaggregate(pwires)

    if (debug)
    {
      plot(pwires, add = T, col = 1:n+1)
      plot(lwires, add = T)
      plot(tlocation, add = T, col = kk)
    }

    # Generate catenary between two consecutive towers
    nlines <- length(pwires)
    Hxy <- vector("list", nlines)
    for (i in 1:nlines)
    {
      # Get the towers for the processing line
      line <- pwires[i,]
      toww <- raster::intersect(tlocation, line)
      tow  <- toww[, c("Z", "deflection")]
      tow$virtual = FALSE

      # Generate virtual towers. Virtual towers are non existing
      # towers that help to prolongate the lines when there is not
      # a second towers
      lin <- raster::intersect(crlwires, line)
      x1 <- sapply(lin@lines, function(x) {
        m <- x@Lines[[1]]@coords
        i <- which.max(m[,1])
        m[i,]
      })
      x1 <- matrix(x1[,which.max(x1[1,])], ncol = 2)

      x2 <- sapply(lin@lines, function(x) {
        m <- x@Lines[[1]]@coords
        i <- which.min(m[,1])
        m[i,]
      })
      x2 <- matrix(x2[, which.min(x2[1,])], ncol = 2)

      vtowers1 <- sp::SpatialPoints(x1, proj4string = proj)
      vtowers2 <- sp::SpatialPoints(x2, proj4string = proj)
      vtowers <- rbind(vtowers1, vtowers2)
      raster::projection(vtowers) <- raster::projection(tow)
      Z <- raster::extract(dtm, vtowers) + mean(tlocation$Z - tlocation$dtm)

      if (anyNA(Z)) stop("Impossible to find DTM value at the edge of the raster. The DTM is not large enought.")

      vtowers <- sp::SpatialPointsDataFrame(vtowers, data.frame(Z, deflection = FALSE, virtual = TRUE))
      tow <- rbind(tow, vtowers)
      #plot(tow, add = T, col = tow$virtual + 1, cex = 2)

      # Order the towers by distance to an arbitrary point so they are in good order in the SPDF
      atower <- toww[1,]
      pmin <- matrix(c(atower@coords[,1] + atower$ux * 5000, atower@coords[,2] + atower$uy * 5000), ncol = 2)
      pmin <- sp::SpatialPoints(pmin, proj4string = proj)
      d <- rgeos::gDistance(pmin, tow, byid = T)
      j <- order(d)
      tow <- tow[j,]
      #plot(tow, add = T, col = "blue")
      #plot(tow, add = T, col = tow$deflection + 2)

      # Remove the virtual towers if not needed.
      # VT are not needed if the they prolongate a line at a deflection point.
      # Special case if there is only a deflection. In this case we are scrapped
      if (sum(tow$deflection) <= sum(!tow$virtual))
      {
        rm = c(FALSE, tow$deflection[-length(tow)] & tow$virtual[-1]) | c(tow$deflection[-1] & tow$virtual[-length(tow)], FALSE)
        tow = tow[!rm,]
      }
      #plot(tiles)
      #plot(tow, add = T, col = "blue")

      # If we have more than a tower we can compute the catenary between two consecutive towers
      # else we do not compute anything
      if (length(tow) > 1)
      {
        hxy <- vector("list", length(tow) - 1)
        for (k in 1:(length(tow) - 1))
        {
          p1 <- list(x = tow@coords[k,1], y = tow@coords[k,2], z = tow$Z[k] - tower.spec$wire.distance.to.top, virtual = tow$virtual[k])
          p2 <- list(x = tow@coords[k + 1, 1], y = tow@coords[k + 1, 2], z = tow$Z[k + 1] - tower.spec$wire.distance.to.top, virtual = tow$virtual[k + 1])
          hxy[[k]] <- catenary(p1, p2, tower.spec$tension)
          hxy[[k]]$virtual <- p1$virtual | p2$virtual
          hxy[[k]]$section <- SECTION
          SECTION <- SECTION + 1
        }

        hxy <- data.table::rbindlist(hxy)
        hxy$ID <- ID
        ID <- ID + 1
        Hxy[[i]] <- hxy
      }
      else
      {
        cc <- list(x = 0, y = 0, z = 0)
        res <- catenary(cc, cc, tower.spec$tension)
        res$virtual = TRUE
        res$section <- 0
        res$ID = 0
        Hxy[[i]] <- res[0,]
      }
    }

    Hxy = data.table::rbindlist(Hxy)
    sp::coordinates(Hxy) = ~x+y
    #rgl::points3d(Hxy@coords[,1],Hxy@coords[,2], Hxy$z)
    raster::projection(Hxy) <- proj
    Hxy <- Hxy[!is.na(sp::over(Hxy, pwlp)),]
    HXY[[kk]] <- Hxy
  }

  HXY <- do.call(rbind, HXY)
  HXY$type = type
  wires <- HXY

  if (debug)
  {
    col <- c("red", "blue", "forestgreen", "darkorchid", "darkorange", "yellow")[wires$ID]
    col[wires$virtual & col == "red"] <- "pink"
    col[wires$virtual & col == "blue"] <- "lightblue"
    col[wires$virtual & col == "forestgreen"] <- "lightgreen"
    col[wires$virtual & col == "darkorchid"] <- "plum"
    col[wires$virtual & col == "darkorange"] <- "goldenrod1"
    col[wires$virtual & col == "yellow"] <- "white"

    plot(raster::extent(dtm), main = paste0("Final extraction"))
    plot(towers, add = T, col = towers$deflection + 1)
    #plot(textent, add = T,  border = textent$deflection + 1)
    plot(wires, col = col, add = T, cex = 0.1)
  }

  # clean wires below ground
  z0 <- dtm[wires]
  invalid <- unique(wires$section[which(wires$z - z0 < 0)])
  wires <- wires[!wires$section %in% invalid,]
  return(wires)
}

# Derivation of Equations for Conductor and Sag Curves of an Overhead Line Based on a Given Catenary Constant
# Alen Hatibovic
# Electrical Engineering and Computer Science 58/1 (2014) 23–27 doi: 10.3311/PPee.6993
catenary = function(p1, p2, c = 1500)
{
  x1 = p1$x
  x2 = p2$x
  y1 = p1$y
  y2 = p2$y
  h1 = p1$z
  h2 = p2$z
  dx = x2-x1
  dy = y2-y1
  dz = h2-h1
  S  = sqrt(dx*dx+dy*dy+dz*dz)
  x  = seq(0, S, by = 1)

  term0 = asinh((h2-h1)/(2*c*sinh(S/(2*c))))
  term1 = x - S/2
  term2 = c*term0
  term3 = S/(2*c)
  term4 = term0

  term5 = sinh(1/(2*c)*(term1+term2))^2
  term6 = sinh(0.5*(term3-term4))^2

  hx = 2*c * (term5 - term6) + h1

  x = seq(x1, x2, length.out = length(hx))
  y = seq(y1, y2, length.out = length(hx))
  z = hx
  return(data.frame(x,y,z))
}

