test_that("get_pmll_estim() worked as anticipated", {
  y_mat <- matrix(1:4, ncol = 2)
  cov_mat <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
  mean_mat <- y_mat - 0.5
  log_fact_vec <- lfactorial(as.vector(y_mat))
  data <- list(y_matrix = y_mat, cov_matrix = cov_mat, mean_matrix = mean_mat,
               log_fact_vec = log_fact_vec)

  set.seed(2)
  expect_equal(round(get_pmll_estim(data), 5), c(-3.50803, -6.68564))
})
