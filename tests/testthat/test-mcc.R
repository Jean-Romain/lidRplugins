context("mcc")

rgdal::set_thin_PROJ6_warnings(TRUE)

file <- system.file("extdata", "Topography.laz", package="lidR")
las = suppressWarnings(readLAS(file, select = "xyzrn"))

test_that("classify_ground mcc works", {

  las <- classify_ground(las, mcc())

  n = names(las@data)

  expect_true("Classification" %in% n)
  expect_equal(unique(las@data$Classification), c(2L, 1L))
  expect_equal(sum(las@data$Classification == 2L), 20777L)
})