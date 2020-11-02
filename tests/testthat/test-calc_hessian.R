test_that("calc_hessian() calculates the hessian correctly", {
  inv_cov <- matrix(c(3, 1, 1, 3), ncol = 2)
  hessian <- round(calc_hessian(c(2, 3), inv_cov, c(-5, -4)), 4)

  expect_equal(hessian, matrix(c(-3.0498, -1.0000, -1.0000, -3.3679), 2))
})
