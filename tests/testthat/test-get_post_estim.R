test_that("get_post_estim() runs as anticipated", {
  dat <- gen_data(sigma = 0.9, phi = 3, seed = 15,
                   patients = 15, time_points_per_patient = 6,
                   alpha_hier = 0.3, beta_hier = 0.2)

  length_gp <- nrow(dat$distance_matrix)
  num_subjects <- length(unique(dat$all_dat$patient))
  data_list <- prepare_data(dat$all_dat, dat$distance_matrix, length_gp,
                            num_subjects, "y", "patient", "time", "position",
                            NULL)

  param_vec <- ifelse(dat$all_true_values %in% c(0.9, 3), 1, 0)
  param_name <- "all"
  subj_idx <- 0
  nsim <- 25

  set.seed(10)
  res <- get_post_estim(data_list, param_vec =  param_vec,
                        param_name =  param_name,
                        subj_index = subj_idx, nsim = nsim)

  res_partial <- get_post_estim(data_list, param_vec =  param_vec - 0.1,
                                param_name =  "alpha",
                                subj_index = 2, nsim = nsim,
                                pmll_vec = res$pmll_vec)

  ## Value was calcuated through previous testing
  expect_equal(round(res$post_value, 3), -1882.646)
  #expect_equal(round(res_partial$post_value, 3), -1860.417)
})
