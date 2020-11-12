test_that("pmmh() works as expected", {
  # Test the base case ------------------------------------------------------

  dat <- gen_data(sigma = 0.9, phi = 3, seed = 16,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.3, beta_hier = 0.2)

  set.seed(10)
  res <- pmmh(dat$all_dat, dat$distance_mat, ndraws = 100, nsim = 3,
              keep_burn_in = TRUE)
  res_df <- evaluate_pmmh(res, dat, burn_in = 51)
  expect_equal(sum(res_df$In_CI), 12)
  expect_equal(sum(res_df$acc_rat), 8.633)
})
