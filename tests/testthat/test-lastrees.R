context("lastrees")

LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
options(lidR.progress = F)
las = lidR::readLAS(LASfile, select = "xyzinr", filter = "-drop_z_below 0 -keep_first -inside 481260 3812921 481300 3812961")

test_that("Hamraz method works", {
  expect_error(lidR::segment_trees(las, hamraz2016()), NA)
})
