test_that("pmmh_single_chain() results match our expectations", {


# Test the base case ------------------------------------------------------

  keep_burn_in <- TRUE
  burn_in <- 51
  outcome <- "y"
  subject_index <- "patient"
  time_index <- "time"
  position_index <- "position"
  ndraws <- 100
  nsim <- 3
  offset_term <- NULL

  dat <- gen_data(sigma = 0.9, phi = 3, seed = 16,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.3, beta_hier = 0.2)

  ## Scalars needed later and prepare list of relevant objects
  length_gp <- nrow(dat$distance_mat)
  num_subjects <- length(unique(dat$all_dat[, subject_index]))
  data_list <- prepare_data(dat$all_dat, dat$distance_mat, length_gp,
                            num_subjects, outcome, subject_index, time_index,
                            position_index, offset_term)
  single_chain_list <- prepare_storage(data_list, ndraws, num_subjects, nsim,
                                       burn_in = burn_in, keep_burn_in = keep_burn_in)

  set.seed(10)
  res <- pmmh_single_chain(single_chain_list)
  res_df <- evaluate_pmmh(res, dat, burn_in = burn_in)
  expect_equal(sum(res_df$In_CI), 11)
  expect_equal(sum(res_df$acc_rat), 7.777)

# Test the addition of a offset term --------------------------------------
  nsim = 10
  offset_term = "offset_vals"

  dat <- gen_data(sigma = 0.5, phi = 2, seed = 14, length_gp = 6,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.4, beta_hier = 0.1,
                  offset_vec = log(seq(1, 3, length.out = 6)))

  ## Scalars needed later and prepare list of relevant objects
  length_gp <- nrow(dat$distance_mat)
  num_subjects <- length(unique(dat$all_dat[, subject_index]))
  data_list <- prepare_data(dat$all_dat, dat$distance_mat, length_gp,
                            num_subjects, outcome, subject_index, time_index,
                            position_index, offset_term)
  single_chain_list <- prepare_storage(data_list, ndraws, num_subjects, nsim,
                                       burn_in, keep_burn_in)

  set.seed(10)
  res <- pmmh_single_chain(single_chain_list)
  res_df <- evaluate_pmmh(res, dat, burn_in = burn_in)
  expect_equal(round(max(res_df$sq_diff), 2), 1.16)
  expect_equal(sum(res_df$In_CI), 10)

  # Test the addition of covariates --------------------------------------
  nsim = 10
  offset_term = NULL
  dat <- gen_data(sigma = 0.5, phi = 2, seed = 14, length_gp = 10,
                  patients = 4, time_points_per_patient = 3,
                  alpha_hier = 0.4, beta_hier = 0.1,
                  num_add_covar = 2)

  ## Scalars needed later and prepare list of relevant objects
  length_gp <- nrow(dat$distance_mat)
  num_subjects <- length(unique(dat$all_dat[, subject_index]))
  data_list <- prepare_data(dat$all_dat, dat$distance_mat, length_gp,
                            num_subjects, outcome, subject_index, time_index,
                            position_index, offset_term)
  single_chain_list <- prepare_storage(data_list, ndraws, num_subjects, nsim,
                                       burn_in, keep_burn_in)

  set.seed(10)
  res <- pmmh_single_chain(single_chain_list)
  res_df <- evaluate_pmmh(res, dat, burn_in = burn_in)
  expect_equal(round(sum(res_df$acc_rat), 2), 8.92)
  expect_equal(sum(res_df$In_CI), 11)
})
