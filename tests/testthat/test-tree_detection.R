context("tree_detection")

LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
las = lidR::readLAS(LASfile, select = "xyz", filter = "-drop_z_below 0")
ctg = lidR::catalog(LASfile)
options(lidR.progress = F)
lidR::opt_progress(ctg) <- FALSE
lidR::opt_chunk_alignment(ctg) <- c(-10,3812970)
lidR::opt_chunk_size(ctg) <- 60
lidR::opt_chunk_buffer(ctg) <- 20

#-------------------------------------------------

test_that("tree_detection mutltichm works with a LAS", {

  ttops = lidR::tree_detection(las, multichm(res = 2, layer_thickness = 2, ws = 5))

  expect_is(ttops, "SpatialPointsDataFrame")
  expect_equal(dim(ttops@data), c(226,2))
  expect_equal(ttops@proj4string, las@proj4string)
})

test_that("tree_detection multichm works with a LAScatalog", {

  ttops = lidR::tree_detection(ctg, multichm(res = 2, layer_thickness = 2, ws = 5))

  expect_is(ttops, "SpatialPointsDataFrame")
  expect_equal(dim(ttops@data), c(226,2))
  expect_equal(ttops@proj4string, ctg@proj4string)
})

#-------------------------------------------------

test_that("tree_detection ptree works with a LAS", {

  ttops = lidR::tree_detection(las, ptrees(k = c(30,15)))

  expect_is(ttops, "SpatialPointsDataFrame")
  expect_equal(dim(ttops@data), c(233,2))
  expect_equal(ttops@proj4string, las@proj4string)
})

test_that("tree_detection ptree works with a LAScatalog", {
  ttops = lidR::tree_detection(ctg, ptrees(k = c(30,15)))

  expect_is(ttops, "SpatialPointsDataFrame")
  expect_equal(dim(ttops@data), c(231,2))
  expect_equal(ttops@proj4string, ctg@proj4string)
})

