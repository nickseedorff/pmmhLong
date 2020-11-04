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
  alpha <- param_vec[data$location_list$alpha]
  beta <- param_vec[data$location_list$beta]
  sigma <- param_vec[data$location_list$sigma]
  phi <- param_vec[data$location_list$phi]
  alpha_hier <- param_vec[data$location_list$alpha_h]
  beta_hier <- param_vec[data$location_list$beta_h]

  ## Commonly used vectors
  subject_index <- data$identifiers$subject_index

  if(param_name %in% c("alpha_h", "beta_h")){
    pmll_est <- pmll_vec
  } else if(param_name %in% c("sigma", "phi", "all")) {

    ## Covariance and precision matrices
    cov_mat <- sigma * exp(-data$distance_mat ^ 2/ phi)
    inv_cov_mat <- qr.solve(cov_mat)

    ## Build mean vector and matrix
    alpha_vec <- alpha[subject_index]
    beta_vec <- beta[subject_index]
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index
    mean_mat <- matrix(mean_vec, nrow = data$length_gp)

    ## Get estimates of the pmll
    pmll_est <- tryCatch(
      get_pmll_estim(data$y_matrix, cov_mat, inv_cov_mat, mean_mat, nsim),
      error  = function(err){
        warning("Error in loop")
        "err"
        })
  } else if(param_name %in% c("alpha", "beta")) {

    ## Covariance and precision matrices
    cov_mat <- sigma * exp(-data$distance_mat ^ 2/ phi)
    inv_cov_mat <- qr.solve(cov_mat)

    ## Build mean vector and matrix
    pat_rows <- which(subject_index == subj_index)
    alpha_vec <- alpha[subject_index[pat_rows]]
    beta_vec <- beta[subject_index[pat_rows]]
    mean_vec <- alpha_vec + beta_vec * data$identifiers$time_index[pat_rows]
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
    pmll_vec <- rep(1e-6, length(data$y_subj_index))
    post_value <- 1e-9
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
