#' Build data matrices for ingestion by gen_pmll_estim()
#'
#' @param data_obj list of objects to make relevant data matrices
#' @param sigma scalar
#' @param phi scalar
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list of matrices
#' @export

gen_data_mat <- function(data_obj, sigma, phi, mean_vec) {
  length_gp <- nrow(data_obj$dist_mat)
  y_matrix <- matrix(data_obj$Y, nrow = length_gp)
  cov_matrix <- sigma * exp(-data_obj$dist_mat ^ 2/ phi)
  mean_matrix <- matrix(mean_vec, nrow = length_gp, ncol = ncol(y_matrix))
  log_fact_vec <- apply(y_matrix, 2, function(x) sum(lfactorial(x)))
  list(y_matrix = y_matrix, cov_matrix = cov_matrix, mean_matrix = mean_matrix,
       log_fact_vec = log_fact_vec)
}

#' Build a dataset for easy testing
#'
#' @param data_obj list of objects to make relevant data matrices
#' @param sigma scalar
#' @param phi scalar
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list of matrices
#' @export

gen_data <- function(sigma, phi, patients, time_points_per_patient,
                     length_gp = 10, seed = 28, alpha_hier, beta_hier) {
  total_obs <- time_points_per_patient * length_gp * patients
  pat_time_combo <- time_points_per_patient * patients
  patient_index <- rep(1:patients, each = time_points_per_patient * length_gp)
  pat_full_index <- rep(1:patients, each = time_points_per_patient)

  ## Alpha should be same for all simulations
  set.seed(29)
  alpha_vals <- rnorm(patients, alpha_hier, sd = 0.5)
  alpha_vec <- alpha_vals[patient_index]

  beta_vals <- rnorm(patients, beta_hier, sd = 0.25)
  beta_vec <- beta_vals[patient_index]

  set.seed(seed)
  time_vec <- rnorm(total_obs)
  ## Distance and covariance matrices
  dist_mat <- as.matrix(dist(0:(length_gp - 1), upper = T, diag = T))
  cov_mat <- sigma * exp(- dist_mat ^ 2/ phi)

  ## Draws from gaussian process
  gauss_process <- t(MASS::mvrnorm(pat_time_combo,
                                   mu = rep(0, length_gp),
                                   Sigma = cov_mat))
  gp_truth <- as.vector(gauss_process)

  ## Linear prediction and generate data
  linear_pred <- alpha_vec + beta_vec * time_vec + gp_truth
  Y <- rpois(total_obs, exp(linear_pred))
  list(Y = Y, dist_mat = dist_mat, gp_truth = gp_truth, time_vec = time_vec,
       alpha_vals = alpha_vals,
       all_true_values = c(alpha_vals, beta_vals, sigma, phi, alpha_hier, beta_hier),
       patients = patients,
       patient_index = patient_index,
       pat_full_index = pat_full_index,
       time_points_per_patient = time_points_per_patient)
}
