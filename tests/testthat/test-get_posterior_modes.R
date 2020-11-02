test_that("get_posterior_modes() runs as expected", {
  y <- matrix(1:4, ncol = 2)
  inv_cov <- qr.solve(matrix(c(1, 0.5, 0.5, 1), ncol = 2))
  mean_mat <- y - 0.5
  mode_res <- get_posterior_modes(y, inv_cov, mean_mat)

  expect_equal(round(mode_res$post_modes, 2),
               matrix(c(-0.38, -0.58, -1.25, -1.79), ncol = 2))
})
