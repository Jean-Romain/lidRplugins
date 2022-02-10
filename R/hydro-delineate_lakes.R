#' Individual waterbodies segmentation
#'
#' Individual surfacic waterbodies segmentation to find the contour of the waterbodies (lakes, rivers).
#' This function is designed to be used on raw point clouds with 1 point per square meter.
#' Point-clouds must be loaded with option \code{filter = "-thin_with_grid 1"}.
#'
#' See the reference for more details. The methods relies on two steps. A raster-based pre-processing
#' step that segments roughtly the waterbodies and a vector-based step that fine accuraly the contours
#' of the lakes. The raster-based step is made of two steps: (a) a very conservative estimation in
#' which the expected output is a  raster that contains all the waterbodies without false positive
#' but the shape of the waterbodies  might be very poorly estimated and (b) an permissive step in which
#' the shapes of the waterbodies are almost correct but with potential false positive. The two steps
#' are then merged internally in  a single accurate raster. Using this raster to speed-up the
#' computation of the vector-based step, the method computes a very accurate contour of the
#' waterbodies.
#'
#' @section Tips for finding adequante parameters:
#' User can apply the function with \code{tol = NULL} to return only the raster-based step output.
#' Also the parameter \code{tol2} can be only 1 numbers to return only one estimation. \code{tol2} is
#' very important and the users must ensure than the conservative step only retains true positive (using
#' a very conservative value such as 1/10000) and the permissive step have a good estimation of the
#' shape even at the cost of many false positives that are anyway removed later using the conservative
#' step.
#' #' \preformatted{
#' # Run only the raster-based conservative step (output is a RasterLayer)
#' r <- delineate_lakes(las, tol = NULL, tol2 = 1/30000)
#' # Run only the raster-based permissive step (output is a RasterLayer)
#' r <- delineate_lakes(las, tol = NULL, tol2 = 2/1000)
#' # Run both raster-based steps (output is a RasterStack)
#' r <- delineate_lakes(las, tol = NULL, tol2 = c(1/30000, 2/1000))
#' }
#'
# @template LAScatalog
#'
#' @section Supported processing options:
#' Supported processing options for a \code{LAScatalog}. For more details see the
#' \link[lidR:LAScatalog-class]{LAScatalog engine documentation}:
#' \itemize{
#' \item \strong{chunk size}: How much data is loaded at once.
#' \item chunk buffer: The buffer size is estimated internally
#' \item \strong{chunk alignment}: Align the processed clusters.
#' \item \strong{progress}: Displays a progress estimate.
#' \item \strong{output_files}: Return the output in R or write each cluster's output in a file. Supported
#' templates are \code{XLEFT}, \code{XRIGHT}, \code{YBOTTOM}, \code{YTOP}, \code{XCENTER}, \code{YCENTER}
#' \code{ID} and, if chunk size is equal to 0 (processing by file), \code{ORIGINALFILENAME}.
#' \item select: 'xyz' are loaded by default.
#' \item filter: "-thin_with_grid 1" is the default.
#' }
#'
#' @param las A LAS or LAScatalog object.
#' @param trim numeric. Bodies with an area smaller than this value are removed. This param reduce
#' oversegmentation of small bodies than may be mistaken with a waterbody
#' @param res numeric. Resolution of the raster use in the pre-processing steps. This parameter is
#' not of major importance
#' @param p numeric. Value between 0 and 1. During the pre-processing raster-based step, cells
#' that have a probability below 'p' to be a lake are not considered in subsequent processing (see
#' references). This parameter is not of major importance
#' @param tol numeric. Positive value. The method tracks the flat regions. This parameter gives the
#' tolerance to the deviation to 'perfectly flat' region. 0 means 0 tolerance so meighboring points
#' must be perfectly aligned exactly on a. A value ranging between 1/10000 and 1/100 might be good.
#' This parameter is of major importance and greatly affect the output.
#' @param tol2 numeric. Positive values. To speed-up the computation a raster-based pre-processing
#' is made. This parameter has the same meaning than \code{tol} but applies to the raster step. It
#' can contains up to 2 values because the preprocessing step allows for two passes (see details and
#' references). This parameter is of major importance and greatly affect the output. A very small
#'  value like 1/1000 is expected.
#' @param th1,th2,k Param of the function \link{shp_hplane} that is used internally. Notice that
#' \code{th3} is set to \code{1-tol}.
#'
#' @return A SpatialPolygons.
#'
#' @references
#' Roussel Jean-Romain, Segmentation of water bodies from ALS data without spectral information (in prep.)
#'
#' @export
#'
#' @examples
#' LASfile <- system.file("extdata", "Topography.laz", package="lidR")
#' las <- readLAS(LASfile, filter = "-thin_with_grid 1")
#'
#' lake <- delineate_lakes(las, trim = 700)
#' plot(las@header)
#' plot(lake, add = TRUE, col = "cornflowerblue")
#'
#' \dontrun{
#' mapview::mapview(lake)
#' }
#' @import methods
#' @import lidR
delineate_lakes <- function(las,  tol = 1/1000, tol2 = c(1/30*tol, 2*tol), trim = 1000, p = 0.5, res = 5, th1 = 25, th2 = 6, k = 10)
{
  UseMethod("delineate_lakes", las)
}

#' @export
delineate_lakes.LAS <- function(las, tol = 1/1000, tol2 = c(1/30*tol, 2*tol), trim = 1000, p = 0.5, res = 5, th1 = 25, th2 = 6, k = 10)
{
  # Defensive programming
  if (!is.null(tol)) lidR:::assert_is_a_number(tol)
  lidR:::assert_is_numeric(tol2)
  lidR:::assert_is_a_number(trim)
  lidR:::assert_is_a_number(p)
  lidR:::assert_is_a_number(res)
  lidR:::assert_is_a_number(th1)
  lidR:::assert_is_a_number(th2)
  lidR:::assert_is_a_number(k)
  lidR:::assert_all_are_in_closed_range(p, 0, 1)
  lidR:::assert_all_are_positive(tol)
  lidR:::assert_all_are_positive(tol2)

  # This is a default empty SpatialPolygons
  empty <- sp::SpatialPolygons(list())
  raster::projection(empty) <- raster::projection(las)

  # Raster estimation
  lidR:::verbose("Precomputing lakes (raster step)...")
  rlakes <- lake_detection_raster(las, tol2, trim, p, res)
  raster::projection(rlakes) <- raster::projection(las)

  # Special case to return the raster and be able to test the input paramters effects
  if (is.null(tol)) return(rlakes)

  # If 1 raster, it is a single pass estimation. If 2 rasters it is a two passes estimation:
  # one conservative estimation with poor contours but few false positives, one tolerant estimation
  # with fair contours but many false positives. Merge the two estimations into a good one
  if (is(rlakes, "RasterLayer")) {
    lidR:::verbose("Computing lakes (vector step)...")
    rlake <- rlakes
  } else if (is(rlakes, "RasterStack")) {
    tmp <- raster::mask(rlakes[[2]], rlakes[[1]])
    ids <- unique(tmp[])
    rlake <- rlakes[[2]]
    rlake[!rlake[] %in% ids] <- NA
  } else {
    stop("Internal error: one or two rasters were expected. Please report this bug.", call. = FALSE)
  }

  # Compute the polygon with the help of the raster
  lidR:::verbose("Computing lakes (vector step)...")
  vlake <- lake_detection_vector(las, rlake, tol, trim, p, th1, th2, k)
  if (is.null(vlake)) return(empty)

  # Simplify the geometry
  vlake <- rgeos::gSimplify(vlake, 1, TRUE)
  vlake@proj4string <- las@proj4string
  return(vlake)
}

#' @export
delineate_lakes.LAScluster <- function(las, tol = 1/1000, tol2 = c(1/30*tol, 2*tol), trim = 1000, p = 0.5, res = 5, th1 = 25, th2 = 6, k = 10)
{
  las <- lidR::readLAS(las)
  if (lidR::is.empty(las)) return(NULL)
  lake <- delineate_lakes(las, tol, tol2, trim, p, res, th1, th2, k)
  return(lake)
}

#' @export
delineate_lakes.LAScatalog <- function(las, tol = 1/1000, tol2 = c(1/30*tol, 2*tol), trim = 1000, p = 0.5, res = 5, th1 = 25, th2 = 6, k = 10)
{
  lidR::opt_select(las) <- "xyz"
  lidR::opt_filter(las) <- "-thin_with_grid 1"
  lidR::opt_chunk_buffer(las) <- ceiling(sqrt(trim))
  options <- list(need_buffer = TRUE)
  output  <- lidR::catalog_apply(las, delineate_lakes, tol = tol, tol2 = tol2, trim = trim, p = p, res = res, th1 = th1, th2 = th2, k = k, .options = options)

  lidR:::verbose("Merging the independent outputs...")
  if (lidR::opt_output_files(las) == "") {
    k <- 1
    n <- length(output)
    for (i in 1:n) {
      m <- length(output[[i]]@polygons)
      if (m == 0) next
      for (j in 1:m) {
        output[[i]]@polygons[[j]]@ID = as.character(k)
        k <- k + 1
      }
    }

    output <- do.call(rbind, output)
    output <- gMerge(output)
    output@proj4string <- las@proj4string
  } else {
    output <- unlist(output)
  }

  return(output)
}

lake_detection_raster <- function(las, tol = 1/1000, trim = 1000, p = 0.5, res = 5)
{
  # Count the number of pass
  n = if (length(tol) == 2L) 2L else 1L
  if (n == 1) lidR:::verbose("Segmentation of the lakes using 1 raster pass")
  if (n == 2) lidR:::verbose("Segmentation of the lakes using 2 raster passes")

  # velox's function `focal` is much faster than the raster's one
  use_velox <- requireNamespace("velox", quietly = TRUE)
  if (!use_velox) message("Install the 'velox' package to speed up this computation")

  # Compute metrics to estimate flatness
  M <- lidR::grid_metrics(las, ~c(stdshapemetrics(X,Y,Z), list(dz = max(Z)-min(Z), sdz = sd(Z), n = length(Z))), res)

  # Smooth with median filter
  if (use_velox) {
    nM <- names(M)
    M  <- velox::velox(M)
    M$medianFocal(3,3, 1:length(nM))
    M  <- M$as.RasterBrick()
    names(M) <- nM
  }

  curvature     <- M$curvature
  anisotropy    <- M$anisotropy
  sphericity    <- M$sphericity
  horizontality <- M$horizontality
  delta         <- M$dz
  density       <- M$n
  sdz           <- M$sdz

  if (!use_velox) {
    ker <- matrix(1, 3, 3)

    curvature     <- raster::focal(curvature,     ker, stats::median, na.rm = T)
    anisotropy    <- raster::focal(anisotropy,    ker, stats::median, na.rm = T)
    sphericity    <- raster::focal(sphericity,    ker, stats::median, na.rm = T)
    horizontality <- raster::focal(horizontality, ker, stats::median, na.rm = T)
    delta         <- raster::focal(delta,         ker, stats::median, na.rm = T)
    density       <- raster::focal(density,       ker, stats::median, na.rm = T)
    sdz           <- raster::focal(sdz,           ker, stats::median, na.rm = T)
  }

  output <- vector("list", n)
  for (i in 1:n) {
    curv <- curvature
    anis <- anisotropy
    sphe <- sphericity
    hori <- horizontality
    delt <- delta
    dens <- density
    sdz. <- sdz

    # Filter pixels that are not flat
    suppressWarnings(curv[curv > tol[i]]     <- NA)
    suppressWarnings(anis[anis < 1 - tol[i]] <- NA)
    suppressWarnings(sphe[sphe > tol[i]]     <- NA)
    suppressWarnings(hori[hori < 1 - tol[i]] <- NA)
    suppressWarnings(delt[delt > 0.1]        <- NA)
    suppressWarnings(dens[is.na(dens)]       <- 0L)
    suppressWarnings(dens[dens > res^2*density(las)*0.4] <- NA)
    suppressWarnings(sdz.[sdz. > 0.05]       <- NA)

    # Compute the proportion of layer that estimated that a pixel is flat
    proba_lake <- is.na(curv) + is.na(delt) + is.na(sdz.) + 4*is.na(hori) + is.na(anis) + is.na(sphe) #+ 2*is.na(hplane)
    proba_lake <- 1 - proba_lake/9
    suppressWarnings(proba_lake[!is.na(dens)] <- 1)

    # Binary layer water body / not water body
    rlake  <- proba_lake >= p
    mlake  <- raster::as.matrix(rlake)

    # Morphology to clean the bodies
    mlake  <- EBImage::closing(mlake,  matrix(1, 3, 3))
    mlake  <- EBImage::fillHull(mlake)

    if (all(mlake)) {
      suppressWarnings(rlake[] <- 1L)
      output[[i]] <- rlake
      next
    }

    # Group the bodies and remove small ones
    mlabel <- EBImage::bwlabel(mlake)
    area   <- table(mlabel)*res^2
    keep   <- which(area >= trim)
    area   <- area[keep]
    keep   <- mlabel %in% as.numeric(names(area))
    mlabel[!keep | mlabel == 0] <- NA

    # Return a RasterLayer
    suppressWarnings(rlake[] <- mlabel)
    output[[i]] <- rlake
  }

  if (n == 1)
    return(output[[1]])
  else
    return(raster::stack(output[[1]], output[[2]]))
}

lake_detection_vector <- function(las, rlake = NULL, tol = 1/1000, trim = 1000, p = 0.5, th1 = 25, th2 = 6, k = 10)
{
  # ---- Pre-process ----
  # If a raster-based estimation of the water
  # bodies is provided we can speed-up the computations

  empty <- NULL

  if (!is.null(rlake)) {
    # Exit early if no water body
    if (all(is.na(rlake[])))
      return(empty)

    # General case: buffer the lake and merge information with point cloud
    mlake <- raster::as.matrix(rlake)
    mlake <- EBImage::dilate(mlake, matrix(1,3,3)) > 0
    mlake <- EBImage::bwlabel(mlake)
    suppressWarnings(rlake[] <- mlake)
    rlake[rlake == 0] <- NA
    las    <- lidR::merge_spatial(las, rlake, "inplake")
    filter <- ~!is.na(inplake)
  } else {
    filter <- NULL
  }

  # ---- Segmentation of the water bodies ----

  # Tree hull is used to compute waterbodies bboxes
  bbox <- lidR::delineate_crowns(las, type = "bbox", attribute = "inplake")
  n    <- length(bbox) + 1
  while (length(bbox) != n) {
    n    <- length(bbox)
    bbox <- sp::disaggregate(rgeos::gUnaryUnion(bbox))
    bbox <- gBboxes(bbox)
  }

  # Remove points classified as water bodies
  lidR:::verbose("Filtering points that belong in flat regions...")
  run <- TRUE
  i   <- 1
  lake <- NULL
  while (run) {
    n1  <- lidR::npoints(las)
    n2  <- sum(lidR:::parse_filter(las, filter))
    las <- lidR::segment_shapes(las, lidR::shp_hplane(th1, th2, 1 - tol, k), "lake", filter = filter)
    las <- lidR::filter_poi(las, !lake)
    n3  <- lidR::npoints(las)
    dn  <- n1 - n3
    dn  <- dn / n2
    run <- dn*100 > 0.1
    lidR:::verbose(glue::glue("Initially {n1} points, processing only {n2}, removed {n1-n3} i.e {round(dn*100, 3)}% of the processed points"))
  }

  lakes <- vector("list", length(bbox))

  # Requiered for the triangulation
  xscale  <- las@header@PHB[["X scale factor"]]
  yscale  <- las@header@PHB[["Y scale factor"]]
  xoffset <- las@header@PHB[["X offset"]]
  yoffset <- las@header@PHB[["Y offset"]]
  scales  <- c(xscale, yscale)
  offsets <- c(xoffset, yoffset)

  for (i in seq_along(bbox)) {
    # Get only the points in the bounding box
    laslake <- lidR::clip_roi(las, bbox[i,])

    # Append points on the edges to avoid some edge artifact
    # (Typically a water body that belongs on an edge of the point cloud)
    X <- append_edgepoints(laslake, raster::extent(bbox[i,]), 2)

    if (is.null(X)) {
      lakes[[i]] <- empty
      next
    }

    # Triangulate the points and trim small triangles
    D <- lidR:::tDelaunay(X, trim = -4, scales = scales, offsets = offsets)

    if (nrow(D) < 4) {
      lakes[[i]] <- empty
      next
    }

    # Compute the contours of the triangulation
    hulls <- tHulls(D, as.matrix(X))

    if (length(hulls) == 0) {
      lakes[[i]] <- empty
      next
    }

    # ---- Post-process and cleaning ----

    # Remove small bodies smaller than 'trim'
    # This removes very dummy  bodies
    #verbose("Trimming small bodies...")
    A      <- raster::area(hulls)
    hulls2 <- hulls[A > trim,]

    if (length(hulls2) == 0) {
      lakes[[i]] <- empty
      next
    }

    if (!suppressWarnings(rgeos::gIsValid(hulls2)))
      hulls2 <- rgeos::gBuffer(hulls2, TRUE, width = 0)

    # Retrieve potential islands inside lakes and create hole
    #verbose("Retrieving islands...")
    hulls2 <- gHole(hulls2)

    # Remove small bodies smaller than 'trim'
    # This removes dummy bodies but keeps the hole above 1000
    #verbose("Trimming small bodies...")
    A      <- raster::area(hulls2)
    hulls2 <- hulls2[A > trim,]

    if (length(hulls2) == 0) {
      lakes[[i]] <- empty
      next
    }

    # Make a SpatialPolygons
    raster::projection(hulls2) <- raster::projection(las)
    lakes[[i]] <- hulls2
  }

  lakes <- Filter(Negate(is.null), lakes)
  if (length(lakes) == 0L) return(NULL)
  if (length(lakes) == 1L) return(lakes[[1]])
  lakes <- do.call(raster::bind, lakes)

  if (rgeos::gIsValid(lakes))
    return(lakes)
  else
    stop("Internal error: invalid geometry.", call. = FALSE)
}

append_edgepoints <- function(las, bbox, space = 2)
{
  xscale  <- las@header@PHB[["X scale factor"]]
  yscale  <- las@header@PHB[["Y scale factor"]]

  X  <- lidR:::coordinates3D(las)
  e  <- bbox
  e@xmax <- e@xmax + lidR:::round_any(0.01, xscale)
  e@xmin <- e@xmin - lidR:::round_any(0.01, xscale)
  e@ymax <- e@ymax + lidR:::round_any(0.01, yscale)
  e@ymin <- e@ymin - lidR:::round_any(0.01, yscale)

  if (e@xmax - e@xmin < space || e@ymax - e@ymin < space)
    return(NULL)

  x  <- lidR:::round_any(seq(e@xmin, e@xmax, space), xscale)
  y  <- lidR:::round_any(seq(e@ymin, e@ymax, space), yscale)
  X2 <- data.frame(X = x,      Y = e@ymin, Z = 0)
  X3 <- data.frame(X = x,      Y = e@ymax, Z = 0)
  X4 <- data.frame(X = e@xmin, Y = y,      Z = 0)
  X5 <- data.frame(X = e@xmax, Y = y,      Z = 0)
  X  <- rbind(X, X2, X3, X4, X5)
  #X  <- data.table::as.data.table(X)
  #X  <- X[!duplicated(X)]
  #X  <- X[!duplicated(X),]
  return(X)
}

gHole <- function(spgeom)
{
  if (length(spgeom) == 1)
    return(spgeom)

  q  <- rgeos::gContains(spgeom, byid = T)
  n  <- ncol(q)
  q  <- q *  !diag(TRUE, n, n)
  rs <- rowSums(q)
  cs <- colSums(q)

  if (all(cs == 0))
    return(spgeom)

  n <- ncol(q)
  polygons <- list()

  for (i in 1:n)
  {
    if (cs[i] == 0 && rs[i] == 0)
    {
      polygons <- append(polygons, spgeom@polygons[i])
    }
    else
    {
      j <- which(q[,i] == 1)

      if (length(j) > 0)
      {
        sp <- rgeos::gDifference(sp::SpatialPolygons(spgeom@polygons[i]), sp::SpatialPolygons(spgeom@polygons[j]))

        if (!is.null(sp))
        {
          sp       <- sp@polygons[[1]]
          sp@ID    <- as.character(i)
          polygons <- append(polygons, sp)
        }
      }
    }
  }

  for (i in seq_along(polygons))
    polygons[[i]]@ID <- as.character(i)

  pp <- sp::SpatialPolygons(polygons)
  return(pp)
}

gBboxes <- function(polygons)
{
  individual_bb <- function(polygon, projection)
  {
    polygon         <- sp::SpatialPolygons(list(polygon), proj4string = projection)
    spatial_bbox    <- as(raster::extent(polygon), "SpatialPolygons")
    spatial_bbox    <- spatial_bbox@polygons[[1]]
    spatial_bbox@ID <- polygon@polygons[[1]]@ID
    return(spatial_bbox)
  }

  polys <- lapply(polygons@polygons, individual_bb, polygons@proj4string)
  spatial_polys <- sp::SpatialPolygons(polys, proj4string = polygons@proj4string)
  return(spatial_polys)
}

gMerge <- function(spgeom)
{
  lidR:::verbose("Computing the cascade union of the polygons...")
  n <- 0
  spgeom = sp::SpatialPolygonsDataFrame(spgeom, data.frame(ID = 1:length(spgeom)))
  while (n != length(spgeom))
  {
    ID <- rgeos::gIntersects(spgeom, byid = TRUE, returnDense = FALSE)
    ID <- unique(ID)
    n  <- length(ID)
    spgeom$ID <- 0L
    for (i in seq_along(ID)) spgeom@data[["ID"]][ID[[i]]] <- i

    spgeom <- stats::aggregate(spgeom, by = list("ID" = spgeom$ID), FUN = function(x) { x[1] })
    spgeom@data[["ID"]] <- NULL
    spgeom@data$Area <- raster::area(spgeom)
  }

  lidR:::verbose("Checking the validity after merging...")
  valid <- suppressWarnings(rgeos::gIsValid(spgeom, TRUE, FALSE))

  if (any(!valid))
  {
    lidR:::verbose(glue::glue("Repairing {sum(!valid)} polygons..."))
    ids <- which(!valid)
    for (id in ids)
    {
      ID <- spgeom@polygons[[id]]@ID
      a  <- sapply(spgeom@polygons[[id]]@Polygons, function(x) x@area)
      h  <- sapply(spgeom@polygons[[id]]@Polygons, function(x) x@hole)
      k  <- !h | (a > 1 & h)
      p  <- sp::Polygons(spgeom@polygons[[id]]@Polygons[k], ID = ID)
      attr(p, "comment") <- paste(as.numeric(h[k]), collapse = " ")
      spgeom@polygons[[id]] <- p
    }
  }

  rgeos::gIsValid(spgeom, TRUE)
  return(as(spgeom, "SpatialPolygons"))
}

tHulls = function(D, X, s = 5)
{
  used <- N <- . <- P1 <- P2 <- NULL

  threads = data.table::getDTthreads()
  data.table::setDTthreads(1L)
  on.exit(data.table::setDTthreads(threads))

  # ==== Get the contours of the partial triangulation ====

  # Order the Delaunay triangulation clockwise
  D = t(apply(D, 1, function(x)
  {
    p1 = c(X[x[1], 1], X[x[1], 2], 1)
    p2 = c(X[x[2], 1], X[x[2], 2], 1)
    p3 = c(X[x[3], 1], X[x[3], 2], 1)
    M = rbind(p1,p2,p3)
    clockwise =  det(M) < 0
    if (clockwise) return(x)
    return(c(x[3], x[2], x[1]))
  }))

  # Compute the vertices (egde of each triangle)
  n  = nrow(D)
  p1 = integer(n*3)
  p2 = integer(n*3)
  for (i in 1:n)
  {
    j = 3*i-2
    tri = as.numeric(D[i,])
    p1[j]   = tri[1]
    p2[j]   = tri[2]
    p1[1 + j] = tri[1]
    p2[1 + j] = tri[3]
    p1[2 + j] = tri[2]
    p2[2 + j] = tri[3]
  }

  vertice = matrix(c(p1, p2), ncol = 2)
  svertice = t(apply(vertice, 1, sort))
  svertice = data.table::as.data.table(svertice)
  vertice = data.table::as.data.table(vertice)
  names(vertice) <- c("p1", "p2")

  # Filter the contour vertices
  notcontour1 = duplicated(svertice)
  notcontour2 = duplicated(svertice, fromLast = TRUE)
  contour     = vertice[!(notcontour1|notcontour2)]
  contour$used = FALSE

  # ==== Extract the polygons from the contour ======

  polygons = list()
  ID = 1L

  # While there are contour there are polygons
  while (nrow(contour) > 0)
  {
    # Initialization
    n       = nrow(contour)
    p1      = numeric(n)
    p2      = numeric(n)

    # The first vertice initiate a polygon
    p1[1]   = contour$p1[1]
    p2[1]   = contour$p2[1]
    data.table::set(contour, 1L, 3L, TRUE)

    # Find the last vertice that close the polygon
    j = which(contour$p1 == p1[1] & contour$used == FALSE)
    s = length(j) == 0L

    if (s) j = which(contour$p2 == p1[1] & contour$used == FALSE)
    if (length(j) == 0L) stop("Internal Error: please report this error (tHull L127)")
    if (length(j) > 1) j = j[1]

    if (s) {
      p1[n] = contour$p1[j]
      p2[n] = contour$p2[j]
    } else {
      p1[n] = contour$p2[j]
      p2[n] = contour$p1[j]
    }

    data.table::set(contour, j, 3L, TRUE)

    # Continue until the we reach the last vertice
    for (i in 2:n)
    {
      j = which(contour$p1 == p2[i - 1] & contour$used == FALSE)
      s = length(j) == 0L

      if (s) j = which(contour$p2 == p2[i - 1] & contour$used == FALSE)
      if (length(j) == 0L) stop("Internal Error: please report this error (tHull L147)")
      if (length(j) > 1) j = j[1]

      if (!s) {
        p1[i] = contour$p1[j]
        p2[i] = contour$p2[j]
      } else {
        p1[i] = contour$p2[j]
        p2[i] = contour$p1[j]
      }

      data.table::set(contour, j, 3L, TRUE)

      if (p2[i] == p1[n]) {
        contour = contour[used == FALSE]
        break
      }
    }

    orderedcontour = data.table::data.table(p1, p2)
    orderedcontour = orderedcontour[orderedcontour$p1 != 0,]

    XX = data.frame(P1 = X[orderedcontour$p1,1], P2 = X[orderedcontour$p1,2])
    data.table::setDT(XX)

    # Fix self intersection
    XX[, N := .N , by = .(P1,P2)]
    XX[c(1, .N), N := 1]
    XX = XX[N == 1][, N := NULL][]

    if (nrow(XX) > 10)
    {
      #cat(ID, " (", nrow(XX), ") ", sep = "")
      p = sp::Polygon(XX)
      p = sp::Polygons(list(p), ID = ID)

      sp = sp::SpatialPolygons(list(p))

      if (!suppressWarnings(rgeos::gIsValid(sp)))
      {
        spsf = sf::st_as_sf(sp)
        spsf = sf::st_make_valid(spsf)
        sp = as(spsf, "Spatial")
        sp = as(sp, "SpatialPolygons")
        p = sp@polygons[[1]]
        p@ID = as.character(ID)
      }

      polygons = append(polygons, p)
      ID = ID + 1
    }
  }

  p = sp::SpatialPolygons(polygons)
  return(p)
}

