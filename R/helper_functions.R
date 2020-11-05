#' Build data matrices for ingestion by gen_pmll_estim()
#'
#' @param distance_mat list of objects to make relevant data matrices
#' @param sigma scalar
#' @param phi scalar
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list of matrices
#' @export

gen_data_mats <- function(distance_mat, sigma, phi, mean_vec) {
  length_gp <- nrow(data_obj$dist_mat)
  cov_matrix <- sigma * exp(-data_obj$dist_mat ^ 2/ phi)
  mean_matrix <- matrix(mean_vec, nrow = length_gp)
  list(cov_matrix = cov_matrix, mean_matrix = mean_matrix)
}

#' Prepare data for loop within pmmh
#'
#' @param data list of objects to make relevant data matrices
#' @param distance_mat matrix, distance between multivariate outcomes
#' @param length_gp scalar, length of the multivariate outcomes
#' @param num_subjects scalar, number of subjects
#'
#' @return list of matrices
#' @export

prepare_data <- function(data, distance_mat, length_gp, num_subjects,
                         outcome, subject_index, time_index,
                         position_index, offset_term) {
  ## Reformatting of inputs
  data$subject_index <- data[, subject_index]
  data$time_index <- data[, time_index]
  data$position_index <- data[, position_index]

  ## Need outcome as a matrix, each column is a multivariate observation
  y_matrix <- matrix(data[, outcome], nrow = length_gp)
  y_subj_index <- matrix(data$subject_index, nrow = length_gp)[1, ]

  ## Remove unneccesary variables, include offset
  vars_to_rm <- c(outcome, subject_index, time_index, position_index)
  non_covariates <- c("outcome", "subject_index", "position_index",
                      "offset_term", "time_index")
  if (!is.null(offset_term)) {
    data$offset_term <- data[, offset_term]
    identifiers <- data[, setdiff(colnames(data), c(offset_term, vars_to_rm))]
  } else {
    identifiers <- data[, setdiff(colnames(data), vars_to_rm)]
  }

  add_x_var <- identifiers[, setdiff(colnames(identifiers), non_covariates)]

  ## Store draws from posterior
  param_names <- c(rep(c("alpha", "beta"), each = num_subjects),
                   c("sigma", "phi", "alpha_h", "beta_h", colnames(add_x_var)))

  ## Data components to pass to get_post_estim
  unique_param <- unique(param_names)
  location_list <- vector("list", length = length(unique_param))
  location_list <- setNames(location_list, unique_param)
  for(i in 1:length(unique_param)) {
    location_list[[i]] <- which(param_names == unique_param[i])
  }

  list(y_matrix = y_matrix, identifiers = identifiers,
       y_subj_index = y_subj_index, add_x_var = add_x_var,
       length_gp = length_gp, distance_mat = distance_mat,
       location_list = location_list, param_names = param_names)
}

#' Build a dataset for easy testing
#'
#' @param data_obj list of objects to make relevant data matrices
#' @param sigma scalar
#' @param phi scalar
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list with data matrix vector or true values
#' @export

gen_data <- function(sigma, phi, patients, time_points_per_patient,
                     length_gp = 10, seed = 28, alpha_hier, beta_hier,
                     offset_vec = NULL, num_add_covar = NULL) {

  set.seed(seed)
  total_obs <- time_points_per_patient * length_gp * patients
  pat_time_combo <- time_points_per_patient * patients
  patient_index <- rep(1:patients, each = time_points_per_patient * length_gp)
  patient_time_index <- rep(1:pat_time_combo, each = length_gp)
  position <- rep(1:length_gp, pat_time_combo)

  ## Draw patient specific slopes and intercepts and time vector
  alpha_vals <- rnorm(patients, alpha_hier, sd = 0.5)
  alpha_vec <- alpha_vals[patient_index]
  beta_vals <- rnorm(patients, beta_hier, sd = 0.25)
  beta_vec <- beta_vals[patient_index]
  time_vec <- rnorm(pat_time_combo)[patient_time_index]

  if (!is.null(num_add_covar)) {
    x_mat <- matrix(rnorm(pat_time_combo * num_add_covar),
                        ncol = num_add_covar)
    x_mat <- x_mat[patient_time_index, ]
    x_beta <- rnorm(num_add_covar)
    colnames(x_mat) <- paste0("add_cov", 1:num_add_covar)
  }

  ## Distance and covariance matrices
  dist_mat <- as.matrix(dist(0:(length_gp - 1), upper = T, diag = T))
  cov_mat <- sigma * exp(- dist_mat ^ 2/ phi)

  ## Draws from gaussian process
  gauss_process <- t(MASS::mvrnorm(pat_time_combo,
                                   mu = rep(0, length_gp),
                                   Sigma = cov_mat))
  gp_truth <- as.vector(gauss_process)
  all_true_values <- c(alpha_vals, beta_vals, sigma, phi, alpha_hier, beta_hier)


  ## Linear prediction and generate data
  if (is.null(offset_vec) & is.null(num_add_covar)) {
    linear_pred <- alpha_vec + beta_vec * time_vec + gp_truth
    y <- rpois(total_obs, exp(linear_pred))
    data <- data.frame(y = y, time = time_vec, position = position,
                       patient = patient_index)
  } else if (is.null(offset_vec) & !is.null(num_add_covar)) {
    linear_pred <- alpha_vec + beta_vec * time_vec + gp_truth +
      x_mat %*% x_beta
    y <- rpois(total_obs, exp(linear_pred))
    data <- data.frame(y = y, time = time_vec, position = position,
                       patient = patient_index, x_mat = x_mat)
    all_true_values <- c(all_true_values, x_beta)

  } else if (!is.null(offset_vec) & is.null(num_add_covar)){
    linear_pred <- alpha_vec + beta_vec * time_vec + gp_truth + offset_vec
    y <- rpois(total_obs, exp(linear_pred))
    data <- data.frame(y = y, time = time_vec, position = position,
                       patient = patient_index, offset_vals = offset_vec)

  } else if (!is.null(offset_vec) & !is.null(num_add_covar)){
    linear_pred <- alpha_vec + beta_vec * time_vec + gp_truth + offset_vec +
      x_mat %*% x_beta
    y <- rpois(total_obs, exp(linear_pred))
    data <- data.frame(y = y, time = time_vec, position = position,
                       patient = patient_index, x_mat = x_mat,
                       offset_vals = offset_vec)
    all_true_values <- c(all_true_values, x_beta)
  }

  list(all_dat = data,
       all_true_values = all_true_values,
       distance_matrix = dist_mat)
}

#' Build a dataset for easy testing
#'
#' @param pmmh_res list of results from the pmmh_function
#' @param data scalar
#' @param burn_in scalar, number of draws to discard
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list with data matrix vector or true values
#' @export

evaluate_pmmh <- function(pmmh_res, gen_data_obj, burn_in = 51) {
  tot_draws <- nrow(pmmh_res$all_draws)
  start <- burn_in + 1
  draws <- pmmh_res$all_draws[start:tot_draws, ]
  accepts <- colMeans(pmmh_res$all_accept[start:tot_draws, ])
  summ_stat <- t(apply(draws, 2, function(x) {
    c(mean(x),
      as.numeric(quantile(x, probs = 0.025)),
      as.numeric(quantile(x, probs = 0.975)),
      sd(x))
    }))

  colnames(summ_stat) <- c("mean", "low_ci", "up_ci", "sd")
  df <- round(data.frame(truth = gen_data_obj$all_true_values,
                         summ_stat, acc_rat = accepts), 3)
  df$In_CI <- as.numeric(df$truth >= df$low_ci &
                           df$truth <= df$up_ci)
  df$sq_diff <- (df$truth - df$mean) ^ 2
  df
}
