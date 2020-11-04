test_that("pmmh() results match our expectations", {


# Test the base case ------------------------------------------------------

  dat <- gen_data(sigma = 0.9, phi = 3, seed = 16,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.3, beta_hier = 0.2)

  set.seed(10)
  res <- pmmh(dat$all_dat, dat$distance_matrix, nsim = 3, ndraws = 100)
  res_df <- evaluate_pmmh(res, dat, burn_in = 51)
  expect_equal(sum(res_df$In_CI), 8)
  expect_equal(sum(res_df$acc_rat), 7.001)

# Test the addition of a offset term --------------------------------------

  dat <- gen_data(sigma = 0.5, phi = 2, seed = 14, length_gp = 6,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.4, beta_hier = 0.1,
                  offset_vec = log(seq(1, 3, length.out = 6)))

  set.seed(10)
  res <- pmmh(dat$all_dat, dat$distance_matrix, nsim = 10, ndraws = 100,
              offset_term = "offset_vals")
  res_df <- evaluate_pmmh(res, dat, burn_in = 51)
  expect_equal(round(max(res_df$sq_diff), 2), 1.14)
  expect_equal(sum(res_df$In_CI), 9)
})
