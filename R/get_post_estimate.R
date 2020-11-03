#' Get values for MH accept step
#'
#' @param gen_data_mat_obj object with list of relevant matrices
#' @param nsim scalar, number of importance samples
#' @return vector
#' @export
#' @examples
#' get_post_estimate(data, param_vec, param_name, param_)

get_post_estimate <- function(data, param_vec, param_name, pat_index,
                              patients, nsim,
                              ll_marg_vec = NULL) {

  ## Parameters to build mean vector
  alpha <- param_vec[1:patients]
  beta <- param_vec[(patients + 1):(2 * patients)]
  sigma <- param_vec[2 * patients + 1]
  phi <- param_vec[2 * patients + 2]
  alpha_hier <- param_vec[2 * patients + 3]
  beta_hier <- param_vec[2 * patients + 4]

  if(param_name %in% c("alpha_h", "beta_h")){
    ll_marg <- ll_marg_vec
  } else if(param_name %in% c("sigma", "phi", "all")) {
    alpha_vec <- alpha[data$patient_index]
    beta_vec <- beta[data$patient_index]
    mean_vec <- alpha_vec + beta_vec * data$time_vec
    data_obj <- gen_data_mat(data, sigma, phi, mean_vec)
    ll_marg <- tryCatch(get_pmll_estim(data_obj, nsim),
                        error  = function(err){
                          warning("Error in loop")
                          "err"
                        })
  } else if(param_name %in% c("alpha", "beta")) {
    pat_rows <- which(data$patient_index == pat_index)
    alpha_vec <- alpha[data$patient_index[pat_rows]]
    beta_vec <- beta[data$patient_index[pat_rows]]
    mean_vec <- alpha_vec + beta_vec * data$time_vec[pat_rows]
    data$Y <- data$Y[pat_rows]
    data_obj <- gen_data_mat(data, sigma, phi, mean_vec)
    ll_marg <- tryCatch(get_pmll_estim(data_obj, nsim),
                        error  = function(err){
                          warning("Error in loop")
                          "err"
                        })
  }

  ## If using all subjects, return full vector
  if (ll_marg[1] == "err") {
    marg_vec <- rep(1e-6, length(data$pat_full_index))
    post_value <- 1e-9
    return(list(marg_vec = marg_vec, post_value = post_value))
  } else if (pat_index == 0) {
    marg_vec <- ll_marg
  } else {
    ll_marg_vec[data$pat_full_index == pat_index] <- ll_marg
    marg_vec <- ll_marg_vec
  }

  post_value <- sum(marg_vec) -
    1.01 * log(sigma) - 0.01 / sigma -
    1.01 * log(phi) - 0.01 / phi -
    sum((alpha - alpha_hier) ^ 2) / 2 / 3 ^ 2 -
    sum((beta - beta_hier) ^ 2) / 2 / 3 ^ 2 -
    alpha_hier ^ 2 / 2 / 3 ^ 2 -
    beta_hier ^ 2 / 2 / 3 ^ 2

  list(marg_vec = marg_vec, post_value = post_value)
}
