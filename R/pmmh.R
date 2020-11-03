#' Get values for MH accept step
#'
#' @param gen_data_obj list of objects
#' @param ndraws number of draws
#' @param nsim number of importance samples
#'
#' @return list
#' @export

pmmh <- function(gen_data_obj, ndraws = 300, nsim = 25) {


  patients <- gen_data_obj$patients

  ## Store draws from posterior
  param_names <- c(rep(c("alpha", "beta"), each = patients), c("sigma", "phi", "alpha_h", "beta_h"))
  pat_idx <- c(1:patients, 1:patients, rep(0, 4))
  num_param <- length(pat_idx)
  accept_values <- param_values <- matrix(NA, nrow = ndraws, ncol = length(param_names))
  ll_res <- vector(length = ndraws)

  ####################################################################
  ## Store distance matrix calculates, matrix determinants, and factorials for future use
  ## Create and object for storing old values
  ## How to check since randomness is different? Can view first
  ##
  ##

  for(i in 2:ndraws){
    if(i == 2) {
      start_time <- Sys.time()
      param_values[1, ] <- ifelse(param_names %in% c("sigma", "phi"), 1, 0)
      ll_obj <- get_post_estim(data = gen_data_obj,
                               param_vec = param_values[1, ],
                               param_name = "all",
                               pat_index = 0,
                               nsim = nsim,
                               patients = patients)

      old_vec <- param_values[1, ]
      pmll_old <- ll_obj$pmll_vec
      ll_old <- ll_obj$post_value
    }

    ## Propose new values
    candidate_vec <- rnorm(num_param, param_values[i - 1, ],
                           c(rep(0.5, patients), 0.1, 0.5, 2, 0.5))
    old_vec <- param_values[i - 1, ]
    accept_vec <- rep(0, num_param)
    for(j in 1:ncol(param_values)) {
      prop_vec <- old_vec
      prop_vec[j] <- candidate_vec[j]

      if(param_names[j] %in% c("sigma", "phi") & prop_vec[j] <= 0) {
        ll_diff <- -100000000
      } else {
        ll_prop_obj <- get_post_estim(data = gen_data_obj,
                                      param_vec = prop_vec,
                                      param_name = param_names[j],
                                      pat_index = pat_idx[j],
                                      nsim = nsim,
                                      patients = patients,
                                      pmll_vec = pmll_old)

        #ll_old <- get_monte_estimate(sigma_old, phi_old)
        ll_diff <- ll_prop_obj$post_value - ll_old
      }

      ## Update beta
      if(log(runif(1)) < ll_diff){
        accept_vec[j] <- 1
        old_vec[j] <- candidate_vec[j]
        pmll_old <- ll_prop_obj$pmll_vec
        ll_old <- ll_prop_obj$post_value
      }
    }
    param_values[i, ] <- old_vec
    accept_values[i, ] <- accept_vec
  }

  ## Create matrix with all draws
  param_names <- c(paste0("alpha", 1:patients), paste0("beta", 1:patients),
                   "sigma", "phi", "alpha_h")

  list(param_names = param_names,
       all_draws = param_values,
       all_accept = accept_values,
       eval_time = difftime(Sys.time(), start_time, units = "mins"))
}
