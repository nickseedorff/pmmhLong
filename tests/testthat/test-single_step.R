test_that("single_step() in newton raphson works as expected", {
  inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
  single_step_res <- round(single_step(c(1, 2), c(2, 3), inv_cov, c(-5, -4)), 4)

  expect_equal(single_step_res, c(3.7402, 4.3120))
})
