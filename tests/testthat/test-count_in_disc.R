X = c(1, 2, 3, 10, 11, 12, 1, 2, 3)
Y = c(3, 2, 1, 5, 7, 6, 8, 9, 7)
Z = c(20, 1, 1, 5, 7, 6, 1, 1, 4)

test_that("count_in_disc works", {

  n <- lidRplugins:::C_count_in_disc(X, Y, c(3, 10), c(2, 5), 2.1, 1L)
  expect_equal(n, c(2,1))
})
