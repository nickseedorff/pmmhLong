test_that("multiplication works", {
  data <- gen_data(sigma = 0.9, phi = 3, seed = 15,
                   patients = 15, time_points_per_patient = 6,
                   alpha_hier = 0.3, beta_hier = 0.2)
  param_vec <- data$all_true_values
  param_name <- "all"
  pat_idx <- 0
  patients <- data$patients
  nsim <- 25

  set.seed(10)
  res <- get_post_estim(data, param_vec, param_name, pat_index = pat_idx,
                    patients, 25)

  ## Value was calcuated through previous testing
  expect_equal(round(res$post_value, 3), -1722.962)
})
