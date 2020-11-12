#' Get values for MH accept step
#'
#' @param dat list of objects
#' @param param_vec vector of proposals
#'
#' @return vector
#' @export

maxim_test <- function(dat, param_vec, num_subjects, nsim = 50) {

  data <- prepare_data(dat$all_dat, dat$distance_matrix, length_gp = 10,
                       num_subjects,
                       outcome = "y", subject_index = "patient",
                       time_index = "time",
                       position_index = "position", offset_term = NULL)


  ## Parameters to build mean vector
  locs <- data$location_list
  alpha <- param_vec[locs$alpha]
  beta <- param_vec[locs$beta]
  sigma <- param_vec[locs$sigma]
  phi <- param_vec[locs$phi]
  alpha_hier <- param_vec[locs$hier_alpha]
  beta_hier <- param_vec[locs$hier_beta]

  ## Commonly used vectors
  subject_index <- data$identifiers$subject_index
  offset_vector <- data$identifiers$offset_term
  add_x_var <- data$add_x_var
  ncol_x_var <- ncol(add_x_var)
  x_var_names <- colnames(add_x_var)
  add_x_var <- as.matrix(add_x_var)
  if (ncol_x_var > 0) {
    x_beta <-param_vec[(locs$hier_beta + 1):length(param_vec)]
  }

  ## Covariance and precision matrices
  cov_mat <- sigma * exp(-data$distance_mat ^ 2/ phi)
  inv_cov_mat <- qr.solve(cov_mat)

  ## Build mean vector and matrix
  alpha_vec <- alpha[subject_index]
  beta_vec <- beta[subject_index]

  ## Add offset if specified
  if (is.null(offset_vector) & ncol_x_var == 0) {
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index
  } else if (!is.null(offset_vector) & ncol_x_var == 0) {
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index +
      offset_vector
  } else if (is.null(offset_vector) & ncol_x_var > 0) {
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index +
      add_x_var %*% x_beta
  } else if (!is.null(offset_vector) & ncol_x_var > 0) {
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index +
      offset_vector + add_x_var %*% x_beta
  }
  mean_mat <- matrix(mean_vec, nrow = data$length_gp)

  ## Get estimates of the pmll
  pmll_vec <- get_pmll_estim(data$y_matrix, cov_mat, inv_cov_mat, mean_mat, nsim)

  post_value <- sum(pmll_vec) -
    #1.01 * log(sigma) - 0.01 / sigma -
    #1.01 * log(phi) - 0.01 / phi -
    sum((alpha - alpha_hier) ^ 2) / 2 / 3 ^ 2 -
    sum((beta - beta_hier) ^ 2) / 2 / 3 ^ 2 -
    alpha_hier ^ 2 / 2 / 3 ^ 2 -
    beta_hier ^ 2 / 2 / 3 ^ 2

  #list(pmll_vec = pmll_vec, post_value = post_value)
  post_value
}
