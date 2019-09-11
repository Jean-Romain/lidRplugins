#' Individual Tree Detection Algorithm
#'
#' This function is made to be used in \link[lidR:tree_detection]{tree_detection}. It implements the
#' LayerStacking algorithm for tree detection based on Ayrey et al (2017) (see references).
#' This function implements only the fisrt part of the method i.e. the detection of the trees.
#'
#' @param start scalar The point cloud is horizontally layered at 1-m intervals starting at 'start'
#' meters above the ground. Default is 0.5 (page 18)
#' @param res scalar Resolution of the CHM computed with a point-to-raster approach. Default is 1 (page 19).
#' @param ws1 scalar Windows radius of the first local maxima use to detected tree tops on the CHM.
#' Default is 3 (page 19)
#' @param buf_size scalar Buffer size placed around each point to build a polygonal buffer around each
#' cluster (figure 1c page 20). Default is 0.5 (page 18)
#' @param ws2 scalar Windows radius of the second local maxima used to detected tree tops on the overlap
#' map. Default is 1.5 (page 20)
#' @param harwood logical. In dense conifer stands with little penetration to the center of the tree,
#' additional weight on the overlap map is given to clusters (page 20). Default is FALSE
#' @param hmin scalar. Point below this threshold cannot initiate a new tree.
#'
#' @references
#' Ayrey, E., Fraver, S., Kershaw, J. A., Kenefic, L. S., Hayes, D., Weiskittel, A. R., & Roth, B. E.
#' (2017). Layer Stacking: A Novel Algorithm for Individual Forest Tree Segmentation from LiDAR Point
#' Clouds. Canadian Journal of Remote Sensing, 43(1), 16–27. https://doi.org/10.1080/07038992.2017.1252907
#' @export
LayerStacking = function(start = 0.5, res = 1, ws1 = 3, ws2 = 1.5, buf_size = 0.5, hardwood = FALSE, hmin = 2)
{
  f = function(las)
  {
    context <- tryCatch({get("lidR.context", envir = parent.frame())}, error = function(e) {return(NULL)})
    lidR:::assert_is_valid_context("tree_detection", "LayerStacking")

    Z <- NULL

    # Page 18: layering ; fig 1a page 20
    las    <- lidR::lasfilter(las, Z > start)
    layers <- round(las@data$Z)
    layers <- split(las@data, layers)

    # Page 19: CHM with resolution of 1
    chm <- lidR::grid_canopy(las, res, lidR::p2r(na.fill = lidR::knnidw(3)))

    # Page 19: Smooth with 3 x 3 m cell
    ksize  <- if (as.integer(3/res) %% 2 == 0) as.integer(3/res) + 1 else as.integer(3/res)
    kernel <- matrix(1, nrow = ksize, ncol = ksize)
    schm   <- raster::focal(chm, w = kernel, fun = function(x){ mean(x, na.rm = TRUE) })

    # Page 19: Local Maximum tree detection
    ttops   <- lidR::tree_detection(schm, lidR::lmf(ws1*2))

    layers_cl    <- LayerStacking_LayerCluster(layers, ttops)         # Page 19: kmeans clustering with ttops used as seed points
    layers_bu    <- LayerStacking_LayerBuffer(layers_cl, buf_size)    # Page 20: Polygonal buffer
    overlap_map  <- LayerStacking_LayerStack(layers_bu, chm, hardwood)

    # Page 20: Local Maximum tree detection
    soverlap_map <- raster::focal(overlap_map, w = kernel, fun = function(x){ mean(x, na.rm = TRUE) })
    ttops        <- lidR::tree_detection(soverlap_map, lidR::lmf(ws2*2, hmin = 0))
    ttops$Z      <- chm[ttops]
    ttops        <- ttops[!is.na(ttops$Z),]
    ttops        <- ttops[ttops$Z >= hmin,]
    return(ttops)
  }

  class(f) <- c("PointCloudBased", "IndividualTreeDetection", "Algorithm", "lidR")
  return(f)
}

LayerStacking_LayerCluster = function(layers, ttops)
{
  for (i in seq_along(layers))
  {
    layer = layers[[i]]

    # Page 18: Density based scanning for the first 3 layer
    # From the original source code since this step is not explained like that
    if (i <= 3)
    {
      if (length(layer) > 1)
      {
        db          <- fpc::dbscan(layer[,1:2], eps = 1.5, MinPts = 5)
        lowlayer    <- layer
        lowlayer$cl <- db$cluster
        lowlayer    <- subset(lowlayer, cl == 0)
        layer       <- lowlayer
      }
    }

    # Remove tree tops lower than the current layer (i is also the height)
    # From the original source code since this step is not explained
    ttops2 <- ttops[ttops$Z > i,]

    # This is expected to be a k-means clustering with the local maxima used as seed points
    # If k-mean does not perform with seed, force to use random seed (not in the paper)
    if (nrow(ttops2) == 0 || nrow(layer) <= nrow(ttops2@coords))
      centers = 1
    else
      centers = ttops2@coords

    group = tryCatch(
    {
      x = stats::kmeans(layer[,1:2], centers)
      x$cluster
    },
    error = function(e)
    {
      warning(paste("kmeans cannot cluster from seed points. Layers", i, "has been removed."), call. = FALSE)
      return(NULL)
    })

    if (!is.null(group))
    {
      layer$cl = group
      layers[[i]] <- layer
    }
  }

  # Remove NULL
  i = which(sapply(layers, function(x){is.null(x$cl)}))

  if (length(i) == 0)
    return(layers)
  else
    return(layers[-i])
}

LayerStacking_LayerBuffer = function(layers, buff)
{
  output = vector(mode = "list", length(layers))

  for (i in 1:length(layers))
  {
    layer = layers[[i]]

    lay_groups1 <- split(layer, layer$cl)

    # To each cluster apply the buffer functions (page 20).
    # It makes a multipart polygon that is the join the 0.5 m radium discs around each point
    buffered_layer1 <- lapply(lay_groups1, LayerBuffer, buf_width = buff)

    # bind the polygons
    buffered_layer1s = do.call(sp::rbind.SpatialPolygons, c(buffered_layer1, makeUniqueIDs = TRUE))

    output[[i]] =  buffered_layer1s
  }

  return(output)
}

LayerStacking_LayerStack = function(layers, layout, hw)
{
  overlap_map = layout
  overlap_map[] <- 0

  for (i in 1:length(layers))
  {
    # What is that ??
    bbb = layers[[i]]
    coords = sapply(bbb@polygons, function(x) sp::coordinates(x@Polygons[[1]]))

    if (length(bbb) > 1 && length(bbb) == length(coords))
      bbb = bbb[!duplicated(coords)]

    layout[] <- 1
    v <- velox::velox(layout)
    spdf <- sp::SpatialPolygonsDataFrame(bbb, data.frame(id = 1:length(bbb)), FALSE)
    v$rasterize(spdf, field = "id", background = NA)
    x <- v$as.RasterLayer()

    layer_raster = x >= 0

    # Adding additional weight to polygons near the top (page 20)
    # This is objectively a weird code
    if (!hw)
    {
      if (i/length(layers) >= .7)
      {
        layer_raster7 = x >= 0
        layer_raster = layer_raster + layer_raster7
      }

      if (i/length(layers) >= .8)
      {
        layer_raster8 = x >= 0
        layer_raster = layer_raster + layer_raster8
      }

      if (i/length(layers) >= .9)
      {
        layer_raster9 = x >= 0
        layer_raster = layer_raster + layer_raster9
      }
    }

    layer_raster[is.na(layer_raster[])] <- 0

    # Construct the overlap map from individual rasters (page 20)
    overlap_map = overlap_map + layer_raster
  }

  return(overlap_map)
}

LayerBuffer <- function(group, buf_width = .6)
{
  pts <- sp::SpatialPoints(group[,1:2])
  buffered_layer = rgeos::gBuffer(pts, width = buf_width)
  buffered_layer = sp::disaggregate(buffered_layer)

  for (i in 1:length(buffered_layer))
  {
    poly = buffered_layer[i]
    outerRings = Filter(function(f){f@ringDir == 1},poly[1]@polygons[[1]]@Polygons)
    outerBounds = sp::SpatialPolygons(list(sp::Polygons(outerRings,ID = 1)))

    if (i == 1)
      buffers = outerBounds
    else
      buffers = sp::rbind.SpatialPolygons(buffers, outerBounds, makeUniqueIDs = TRUE)
  }

  return(buffers)
}