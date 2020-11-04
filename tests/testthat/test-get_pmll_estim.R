test_that("get_pmll_estim() worked as anticipated", {
  y_mat <- matrix(1:4, ncol = 2)
  cov_mat <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
  inv_cov_mat <- qr.solve(cov_mat)
  mean_mat <- y_mat - 0.5

  set.seed(2)
  obs_val <- round(get_pmll_estim(y_mat, cov_mat, inv_cov_mat,
                                  mean_mat, nsim = 500), 5)
  expect_equal(obs_val, c(-3.50803, -6.68564))
})
