#' Get values for MH accept step
#'
#' @param data list of objects
#' @param param_vec vector of proposals
#' @param param_name char, name of parameter proposed
#' @param subj_index scalar, index of subject
#' @param nsim scalar, number of importance samples
#' @param pmll_vec vector of previous pmll estimates
#'
#' @return vector
#' @export

get_post_estim <- function(data, param_vec, param_name, subj_index, nsim,
                           pmll_vec = NULL) {

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


  if(param_name %in% c("hier_alpha", "hier_beta")){
    pmll_est <- pmll_vec
  } else if(param_name %in% c("sigma", "phi", "all", x_var_names)) {

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
    pmll_est <- tryCatch(
      get_pmll_estim(data$y_matrix, cov_mat, inv_cov_mat, mean_mat, nsim),
      error  = function(err){
        warning("Error in loop")
        "err"
        })
  } else if (param_name %in% c("alpha", "beta")) {

    ## Covariance and precision matrices
    cov_mat <- sigma * exp(-data$distance_mat ^ 2/ phi)
    inv_cov_mat <- qr.solve(cov_mat)

    ## Build mean vector and matrix
    pat_rows <- which(subject_index == subj_index)
    alpha_vec <- alpha[subject_index[pat_rows]]
    beta_vec <- beta[subject_index[pat_rows]]

    ## Add offset if specified
    if (is.null(offset_vector) & ncol_x_var == 0) {
      mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index[pat_rows]
    } else if (!is.null(offset_vector) & ncol_x_var == 0) {
      mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index[pat_rows] +
        offset_vector[pat_rows]
    } else if (is.null(offset_vector) & ncol_x_var > 0) {
      mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index[pat_rows] +
        add_x_var[pat_rows, ] %*% x_beta
    } else if (!is.null(offset_vector) & ncol_x_var > 0) {
      mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index +
        offset_vector[pat_rows] + add_x_var[pat_rows, ] %*% x_beta
    }

    mean_mat <- matrix(mean_vec, nrow = data$length_gp)

    ## Get estimates of the pmll for only individuals where it will change
    pmll_est <- tryCatch(
      get_pmll_estim(data$y_matrix[, data$y_subj_index == subj_index], cov_mat,
                     inv_cov_mat, mean_mat, nsim),
      error  = function(err){
        warning("Error in loop")
        "err"
        })
  }

  ## If using all subjects, return full vector
  if (pmll_est[1] == "err") {
    pmll_vec <- rep(-1e12, length(data$y_subj_index))
    post_value <- -1e12
    return(list(pmll_vec = pmll_vec, post_value = post_value))
  } else if (subj_index == 0) {
    pmll_vec <- pmll_est
  } else {
    pmll_vec[data$y_subj_index == subj_index] <- pmll_est
  }

  post_value <- sum(pmll_vec) -
    1.01 * log(sigma) - 0.01 / sigma -
    1.01 * log(phi) - 0.01 / phi -
    sum((alpha - alpha_hier) ^ 2) / 2 / 3 ^ 2 -
    sum((beta - beta_hier) ^ 2) / 2 / 3 ^ 2 -
    alpha_hier ^ 2 / 2 / 3 ^ 2 -
    beta_hier ^ 2 / 2 / 3 ^ 2

  list(pmll_vec = pmll_vec, post_value = post_value)
}
