test_that("calc_gradient() calculates the correct gradient value", {
  inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
  gradient <- round(calc_gradient(c(1, 2), c(2, 3), inv_cov, c(-5, -4)), 4)

  expect_equal(gradient, c(0.5752, 0.7571))
})
