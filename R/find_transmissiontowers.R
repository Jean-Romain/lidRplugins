#' Individual transmission tower detection
#'
#' Individual transmission tower detection function that find the positions of the transmission
#' towers. The method is supervised by a map of the electric network and the tower types.
#'
#' @param las An object of class LAS with absolute elevations or a LAScatalog.
#' @param powerline A \code{SpatialLines*} that map the electrical network accurately
#' @param type character. One of "waist-type", "waist-type-small", "double-circuit" according to
#' \href{http://www.hydroquebec.com/learning/transport/types-pylones.html}{Hydro-Quebec}. Can also
#' be a list with custom specifications. See \link{get_tower_spec}.
#' @param buffer numeric. The \code{SpatialLines*} will be buffered internally to catch the powerlines
#' and the transmission towers. The buffer must ensure to catch all the powerlines.
#' @param dtm \code{RasterLayer}. Because the algorithm relies on absolute elevation a DTM is
#' requirered to compute the relative elevations.
#' @param debug logical. Plot the different steps of the algorithm so one can try to figure out what
#' is going wrong.
#'
#' @return A \code{SpatialPointDataFrame} with several attributes. \code{Z} the elevation of the tower,
#' \code{dtm} the elevation of the bottom of the tower aligned with the top, \code{theta} the angle
#' of the tower with the x axis in radian, \code{ux, uy} the directional vectors, \code{deflection}
#' tells if a given tower is on a deflection (deflection towers are found twice by design) and
#' \code{type} the type name.
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
#' dtm <- grid_terrain(las, 2, tin())
#'
#' towers <- find_transmissiontowers(las, network, dtm, "waist-type")
#'
#' plot(las@header)
#' plot(towers, add = T, col = towers$deflection + 1)
#' arrows(
#'    towers@coords[,1],
#'    towers@coords[,2],
#'    towers@coords[,1] + 100 * towers$ux,
#'    towers@coords[,2] + 100 * towers$uy,
#'    length = 0.05,
#'    col = towers$deflection + 1)
#'
#' plot(las) %>% add_treetops3d(towers, radius = 5)
#' }
#' @family electric network
#' @export
find_transmissiontowers = function(las, powerline, dtm, type = c("waist-type", "double-circuit"), buffer = 125, debug = FALSE)
{
  UseMethod("find_transmissiontowers", las)
}

#' @export
find_transmissiontowers.LAS = function(las, powerline, dtm, type = c("waist-type", "double-circuit"), buffer = 125, debug = FALSE)
{
  lidR:::assert_is_all_of(powerline, "SpatialLinesDataFrame")
  lidR:::assert_is_all_of(dtm, "RasterLayer")
  lidR:::assert_all_are_positive(buffer)

  if (debug)
  {
    opar = graphics::par("mfrow")
    graphics::par(mfrow = c(2,3))
    on.exit(graphics::par(mfrow = opar))
  }

  # The name 'las' will be used later. Keep original object untouched
  olas <- las

  # Get the spec of the transmission towers
  tower.spec <- get_tower_spec(type)

  # Crop the lines to the extent of the las
  pwll <- raster::crop(powerline, raster::extent(las))

  if (debug)
  {
    plot(las@header, main = paste0("Raw powerline network"))
    plot(pwll, add = T, col = 1:length(pwll))
  }

  pwll <- gJoinLines(pwll, 2)
  pwll <- rgeos::gSimplify(pwll, 40)

  # Split each segment/section of the powerline
  # This allow to support powerline deflection
  spwll <- gSplitLines(pwll)
  spwll <- gElongateLines(spwll, 25)
  if (length(spwll) != length(spwll))
    stop("Internal error: different sizes for spatial objects", call. = TRUE)

  # Transform each segment/section into pol ygon
  spwlp <- rgeos::gBuffer(spwll, width = buffer, byid = T, capStyle = 'SQUARE')

  if (debug)
  {
    plot(las@header, main = "Post-processed lines and buffers")
    #plot(as(spwll, "SpatialPoints"), add = T)
    plot(spwll, add = T, col =  1:length(spwll))
    plot(spwlp, add = T, border =  1:length(spwlp), lty = 3)
  }

  if (debug)
  {
    plot(olas@header, main = paste0("Tower candidates and corrected candidates"))
  }

  # Loop on each segment
  output <- vector("list", length(spwll))
  for (k in 1:length(output))
  {
    # Keep only the section k
    pwll <- spwll[k,]
    pwlp <- spwlp[k,]

    # Clip the las to work only within the buffer of the section of the powerline
    # TODO: we need only XYZ, we can save memory with a better clipping
    las <- lidR::lasclip(olas, pwlp)
    if (!is(las, "LAS")) stop("Internal error: object is not a LAS.")
    if (lidR::is.empty(las)) stop("Internal error: object LAS is empty.")

    # Merge the DTM, it gonna be useful later on
    las <- lidR::lasmergespatial(las, dtm, "dtm")

    # Compute the orientation of the wires and towers
    orientation <- sp::coordinates(pwll)[[1]][[1]]
    orientation <- lidR:::fast_eigen_values(orientation)$coef
    angle <- atan(orientation[2,1]/orientation[1,1])

    # Find candidate location at being a transmission tower
    towers <- tower.candidates(las, dtm, tower.spec, angle)

    if (debug)
    {
      plot(towers, add = T, col = "gray40")
      text(towers@coords[,1], towers@coords[,2]+40, 1:length(towers), cex = 0.8, col = k)
      #plot(las) %>% add_treetops3d(towers, radius = 5)
    }

    # Add informations in table of attribute about tower orientation
    towers$theta <- round(angle,3)
    towers$ux <- round(orientation[1,1],3)
    towers$uy <- round(orientation[2,1],3)

    # The canditate towers are innaccurate with false positive. The following
    # steps aims to clean that by rectifying the positionning of the towers that are no centered
    # on the towers (the ears of the towers are actually detected) and we clear false positives
    rtowers <- tower.rectification(las, towers, tower.spec, angle, dtm)

    if (debug)
    {
      plot(rtowers, add = T, col = k)
      graphics::text(rtowers@coords[,1], rtowers@coords[,2]+40, 1:length(rtowers), cex = 0.8, col = k)
      #plot(las) %>% add_treetops3d(rtowers, radius = 5)
    }

    # Clean some remaining false positive in deflection
    # (I don't remember which case it covers)
    if (length(output) > 1) {
      pwlp2 <- rgeos::gBuffer(pwlp, width = -10, byid = T, capStyle = 'SQUARE')
      keep = rgeos::gWithin(rtowers, pwlp2, byid = T)
      towers <- rtowers[as.logical(keep),]
    } else {
      towers <- rtowers
    }

    output[[k]] <- towers
  }

  # Merge the different segments into a single object
  ptowers <- lapply(output, function(x) x[[1]])
  ntowers <- sum(sapply(ptowers, length))

  if (ntowers > 0)
  {
    ptowers <- do.call(rbind, output)

    points_matrix <- rgeos::gWithinDistance(ptowers, dist = tower.spec$length[2]/2, byid = TRUE)
    diag(points_matrix) <- NA
    v <- colSums(points_matrix, na.rm = TRUE) == 0
    ptowers$deflection = !v

    any_deflection = any(!v)

    if (debug & any_deflection)
    {
      plot(olas@header, main = "All towers before deflection correction")
      plot(ptowers, add = T, col = ptowers$deflection + 1)
      #plot(textent, add = T,  border = tlocation$deflection + 1)
      graphics::arrows(ptowers@coords[,1], ptowers@coords[,2], ptowers@coords[,1] + 100 * ptowers$ux,  ptowers@coords[,2] + 100 * ptowers$uy, length = 0.05, col = ptowers$deflection + 1)
    }

    # last correction for deflection
    in_several_lines = rgeos::gContains(spwlp, ptowers, byid = T)
    remove = rep(FALSE, length(ptowers))

    for (k in 1:length(output))
    {
      pwll <- spwll[k,]

      orientation <- sp::coordinates(pwll)[[1]][[1]]
      orientation <- lidR:::fast_eigen_values(orientation)$coef
      angle <- atan(orientation[2,1]/orientation[1,1])

      in_this_lines = in_several_lines[,k]

      good_angle = abs(ptowers$theta - angle) < 10e-3
      remove[in_this_lines & !good_angle & !ptowers$deflection] <- T
    }

    ptowers = ptowers[!remove,]
    ptowers$type = tower.spec$name

    if (debug)
    {
      plot(olas@header, main = "Final towers")
      plot(ptowers, add = T, col = ptowers$deflection + 1)
      #plot(textent, add = T,  border = tlocation$deflection + 1)
      graphics::arrows(ptowers@coords[,1], ptowers@coords[,2], ptowers@coords[,1] + 100 * ptowers$ux,  ptowers@coords[,2] + 100 * ptowers$uy, length = 0.05, col = ptowers$deflection + 1)
    }

    return(ptowers)
  }
  else
  {
    return(output[[1]])
  }
}

#' @export
find_transmissiontowers.LAScluster = function(las, powerline, dtm, type = c("waist-type", "double-circuit"), buffer = 125, debug = FALSE)
{
  bbox <- raster::extent(las)
  las <- lidR::readLAS(las)
  if (lidR::is.empty(las)) return(NULL)

  # pos and extent enforced to TRUE to guarantee to remove buffer properly
  output <- find_transmissiontowers(las, powerline, dtm, type, buffer)
  output <- raster::crop(output, bbox)
  return(output)
}

#' @export
find_transmissiontowers.LAScatalog = function(las, powerline, dtm, type = c("waist-type", "double-circuit"), buffer = 125, debug = FALSE)
{
  pwrlp <- rgeos::gBuffer(powerline, width = buffer)
  ctg <- lidR::catalog_intersect(las, pwrlp)
  las$processed <- FALSE
  las$processed[row.names(las) %in% row.names(ctg)] <- TRUE

  options = list(need_buffer = TRUE)
  output <- lidR::catalog_sapply(las, find_transmissiontowers, powerline = powerline, type = type, buffer = buffer, dtm = dtm, .options = options)
  return(output)
}

# Using a point cloud, a DTM, the tower specification and knowledge about the orientation
# of the tower, applies a local maximum filter to find which point are likely to be a tower.
# We use the normalized + raw data + soothed raw data because the LMF is prone at missing some
# towers. Using different post process allows to find more towers (i.e. more false positive) but
# the will be cleaned later. What matter is not having false negative.
tower.candidates = function(las, dtm, tower.spec, angle)
{
  splas = lasfiltersurfacepoints(las, 1)

  # We work with relative elevations both in raw + smoothed data
  #nlas <- lidR::lasnormalize(las, dtm)
  ssplas <- lidR::lassmooth(splas, 10)

  # Keep the top 5 m below the towers. With this we are sure to remove most of the noise
  Z <- NULL
  sub <-  lidR::lasfilter(splas, Z > dtm + tower.spec$height[1] - 5)
  #nsub <- lidR::lasfilter(nlas, Z > tower.spec$height[1] - 5)
  ssub <- lidR::lasfilter(ssplas, Z > dtm + tower.spec$height[1]/2)

  # Find the local max using an oriented windows using both raw an normalized data
  # This because in both we could miss some some towers but not the same.
  rtowers <- lidR::find_localmaxima(sub, c(200, tower.spec$length[2]*1.2, angle))
  #ntowers <- lidR::local_maximum(nsub, c(150, tower.spec$length[2]*1.2, angle))
  stowers <- lidR::find_localmaxima(ssub, c(200, tower.spec$length[2]*1.2, angle))

  #ntowers$Z <- ntowers$Zref
  #ntowers$Zref <- NULL
  stowers$Z <- stowers$Zraw
  stowers$Zraw <- NULL
  #plot(las) %>% add_treetops3d(rtowers, radius = 7)
  #plot(las) %>% add_treetops3d(ntowers, radius = 7)
  #plot(las) %>% add_treetops3d(stowers, radius = 4)

  if (length(rtowers) == 0 && length(stowers) == 0)
    return(rtowers)

  # Keep only one tower if duplicates
  towers <- rbind(stowers, rtowers)
  towers <- towers[!duplicated(towers@data),]
  #plot(las) %>% add_treetops3d(towers, radius = 7)

  # We keep only one tower if a tower has been found twice at two close locations
  points_matrix <- rgeos::gWithinDistance(towers, dist = tower.spec$length[2]*1.2, byid = TRUE)
  points_matrix[lower.tri(points_matrix, diag = TRUE)] <- FALSE
  v <- rowSums(points_matrix) == 0
  towers = towers[v, ]
  #plot(las) %>% add_treetops3d(towers, radius = 7)

  # Get some loose on the bottom because some might be underestimated
  rm = towers$Z - towers$dtm >= tower.spec$height[1] - 2 & towers$Z - towers$dtm <= tower.spec$height[2]

  towers = towers[rm, ]
  #plot(las) %>% add_treetops3d(towers, radius = 7)

  #towers = raster::crop(towers, raster::extent(towers) - 2)

  return(towers)

  #plot(towers)
  #text(towers@coords[,1], towers@coords[,2]+25, 1:length(towers), cex = 0.8)
  #plot(las) %>% add_treetops3d(towers, radius = 7, col = colSums(points_matrix)+1)
  #plot(las) %>% add_treetops3d(towers, radius = 7, col =  colSums(points_matrix)+1)
}

tower.rectification <- function(las, towers, tower.spec, angle, dtm)
{
  Z <- NULL

  buffer.towers <- rgeos::gBuffer(towers, width = tower.spec$length[2]/2*1.3, byid = T)
  las2   <- lidR::lasfilter(las, Z > dtm + 2)

  coords <- vector("list", length(buffer.towers))
  for (i in 1:length(buffer.towers))
  {
    Zbottom <- raster::extract(dtm, towers[i,])
    sub2 <- lidR::lasclip(las2, buffer.towers[i,])
    sub2 <- lidR::lasfilter(sub2, Z > Zbottom)
    coords[[i]] <- tower.correction(sub2, angle, tower.spec, Zbottom)
  }

  coords <- data.table::rbindlist(coords)
  rm     <- coords$Tower

  # No towers, return an empty SpatialPolygonsDataFrame and jump to next iteration
  if (all(!rm))
  {
    out <- towers[0,]
    out@bbox <- las@bbox
    return(out)
  }

  # Finalize the positionning of the towers in a SpatialPointsDataFrame
  coordinates <- as.matrix(coords[rm, 1:2])
  rectified.towers <- towers[rm,]
  rectified.towers@coords <- coordinates
  rectified.towers$Z = coords$Z[rm]
  return(rectified.towers)
}

tower.correction <-  function(las, angle, tower.spec, Zbottom)
{
  dtm <- Z <- NULL

  Xm <- mean(las$X[las$Z > Zbottom + 2])
  Ym <- mean(las$Y[las$Z > Zbottom + 2])
  Zm <- max(las$Z)

  z <- las$Z
  nz <- las$Z - las$dtm

  # Test if the distibution is almost continuous on Z
  # A tower is a continuous structure from the ground to the top
  h <- graphics::hist(z, breaks = seq(min(floor(Zbottom) + 5, floor(min(z))), ceiling(max(z)), 1), plot = FALSE)
  G <- sum(h$counts == 0L) <= 10

  # Test if the distibution does not have too much point on the bottom
  # A tower is a continuous structure from the ground to the top
  h <- graphics::hist(nz[nz > 5], breaks = seq(min(5, floor(min(nz))), ceiling(max(nz)), 1), plot = FALSE)
  J <- cumsum(h$density)[as.integer(length(h$density)/2)] < 0.53

  # Test if the distibution is almost continuous on X after reorientation
  # A tower is a continuous structure horizontally
  a <- angle + pi/2
  rot <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol = 2)
  coords <- as.matrix(lidR:::coordinates(lasfilter(las, Z > Zm - tower.spec$wire.distance.to.top - 2)))
  zero <- sp::bbox(las)[,1]
  coords[,1] <- coords[,1] - zero[1]
  coords[,2] <- coords[,2] - zero[2]
  coords <- coords %*% rot
  X <- lidR:::round_any(coords[,1], las@header@PHB[["X scale factor"]])
  #Y <- lidR:::round_any(coords[,2], las@header@PHB[["Y scale factor"]])
  #las@data[["X"]] <- X
  #las@data[["Y"]] <- Y
  #las <- lidR:::lasupdateheader(las)
  fminX <- floor(min(X))
  cmaxX <- ceiling(max(X))
  h <- graphics::hist(X, breaks = seq(fminX, max(fminX + tower.spec$length[1], cmaxX), 1), plot = FALSE)
  K <- sum(h$counts == 0L) <= tower.spec$length[1]/2

  # Test if the area covered in small
  A <- area(las) <= (pi * (tower.spec$length[2]/2*1.3)^2)*0.8

  # There is at most 1/4 test that says it's not a tower: its a tower
  is.tower = J+G+2*K+A >= 4

  S = Zm - Zbottom >= tower.spec$height[1]
  if (!S) is.tower = FALSE

  ret <- list(X = Xm, Y = Ym, Z = Zm, Tower = is.tower)
  return(ret)
}

#' Get specification of a tower
#'
#' Return the specifications of a given tower type i.e. size range, width range, tension, number of
#' wire. If \code{type} is contains the specifications of a non supported type it checks the validity
#' of the specification
#'
#' @param type one of \code{"waist-type"} \code{"double-circuit"} or \code{"waist-type-small"} or a
#' \code{list}
#'
#' @return  A list
#'
#' @examples
#' specs = get_tower_spec("waist-type")
#'
#' # Create new specs
#' new_specs = specs
#' new_specs$width = c(20, 25)
#'
#' # Validate this specs
#'  get_tower_spec(new_specs)
#' @export
#' @family electrical network
get_tower_spec = function(type)
{
  waist.type = list(
    name = "waist-type",
    length = c(38,40),
    width = c(15, 20),
    height = c(32,64),
    wires = 3L,
    wire.layers = 1L,
    wire.distance = 0,
    wire.distance.to.top = 10,
    tension = 1300)

  double.circuit = list(
    name = "double-circuit",
    length = c(15,20),
    width = c(10, 13),
    height = c(40,62),
    wires = 2L,
    wire.layers = 3L,
    wire.distance = 7 ,
    wire.distance.to.top = 8,
    tension = 1400)

  waist.type.small = list(
    name = "waist-type-small",
    length = c(16,22),
    width = c(7, 10),
    height = c(28,40),
    wires = 3L,
    wire.layers = 1L,
    wire.distance = 0,
    wire.distance.to.top = 5,
    tension = 1800)

  if (!is.list(type))
  {
    type <- match.arg(type, c("waist-type", "waist-type-small", "double-circuit"))

    if (type == "waist-type")
      tower.spec <- waist.type
    else if (type == "waist-type-small")
      tower.spec <- waist.type.small
    else if (type == "double-circuit")
      tower.spec <- double.circuit
    else
      stop("This type of transmission tower does not exist", call. = FALSE)

    return(tower.spec)
  }
  else
  {
    if (!all.equal(names(type), names(waist.type)))
      stop("Invalid definition of a tower type: incorrect element names names")

    typeref = sapply(waist.type, typeof)
    typedat = sapply(type, typeof)

    if (!all.equal(typeref, typedat))
      stop("Invalid definition of a tower type: incorrect element types")

    lengthref = sapply(waist.type, length)
    lenthdat = sapply(type, length)

    if (!all.equal(lengthref, lenthdat))
      stop("Invalid definition of a tower type: incorrect element sizes")

    return(invisible(type))
  }
}

# Derivation of Equations for Conductor and Sag Curves of an Overhead Line Based on a Given Catenary Constant
# Alen Hatibovic
# Electrical Engineering and Computer Science 58/1 (2014) 23–27 doi: 10.3311/PPee.6993
catenary = function(p1, p2, c = 1500) {
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

gSplitLines <- function(sl)
{
  ccs <- sp::coordinates(sl)

  out = vector("list", length(ccs))
  for (j in 1:length(ccs))
  {
    cc = ccs[[j]]
    if (length(cc) > 1) stop("Internal error length(cc) > 1 in gSplit")
    cc = cc[[1]]

    outputlist <- vector("list", nrow(cc) - 1)
    i <- 1
    while (i < (nrow(cc)))
    {
      coords1 <- cc[i,]
      coords2 <- cc[i+1,]
      bind <- rbind(coords1, coords2)
      outputlist[[i]] <- sp::Lines(list(sp::Line(bind)), as.character(i))
      i <- i+1
    }

    out[[j]] <- sp::SpatialLines(outputlist, proj4string = sl@proj4string)
  }

  out <- do.call(rbind, out)
  return(out)
}


gElongateLines <- function(sl, l = 25)
{
  ccs <- sp::coordinates(sl)
  ccs2 = ccs
  for (j in 1:length(ccs))
  {
    cc = ccs[[j]]
    if (length(cc) > 1) stop("Internal error: length(cc) > 1")
    cc = cc[[1]]
    if (nrow(cc) > 2) stop("Internal error: spatial lines are not made of simple segments")

    orientation <- lidR:::fast_eigen_values(cc)$coef
    angle <- atan(orientation[2,1]/orientation[1,1])

    line1 = sp::Line(cc)

    xs = if (cc[1,1] < cc[2,1]) -1 else 1
    if (sign(angle) > 0)
      ys = if (cc[1,2] < cc[2,2]) -1 else 1
    else
      ys = if (cc[1,2] < cc[2,2]) 1 else -1

    cc[1,1] = cc[1,1] + xs * l * cos(angle)
    cc[2,1] = cc[2,1] - xs * l * cos(angle)
    cc[1,2] = cc[1,2] + ys * l * sin(angle)
    cc[2,2] = cc[2,2] - ys * l * sin(angle)

    line2 = sp::Line(cc)

    if (sp::LineLength(line2) <= sp::LineLength(line1))
      stop("Internal error: line not elongated in gElongate")

    ccs2[[j]] <- sp::Lines(line2, ID = as.character(j))
  }

  out = sp::SpatialLines(ccs2, proj4string = sl@proj4string)
  return(out)
}

gJoinLines = function(sl, th = 2)
{
  cc <- sp::coordinates(sl)
  cc <- lapply(cc, function(x) { do.call(rbind, x) })
  cc <- do.call(rbind, cc)
  sp <- sp::SpatialPoints(cc)
  m <- rgeos::gWithinDistance(sp, dist = 5, byid = T)
  m[upper.tri(m, diag = TRUE)] <- FALSE
  join <- which(m, arr.ind = T)

  if (nrow(join) > 1) stop("Internal error: to many lines to join")
  if (nrow(join) == 0) return(sl)

  cc <- cc[as.numeric(join),]
  xm <- mean(cc[,1])
  ym <- mean(cc[,2])

  for (i in 1:length(sl))
  {
    l <- sl@lines[[i]]@Lines[[1]]@coords
    u <- l[,1] %in% cc[,1] & l[,2] %in% cc[,2]
    sl@lines[[i]]@Lines[[1]]@coords[u,1] <- xm
    sl@lines[[i]]@Lines[[1]]@coords[u,2] <- ym
  }

  sl2 <- rgeos::gLineMerge(sl)
  sl2 <- sp::disaggregate(sl2)
  return(sl2)
}
