#' Run pmmh for a single chain
#'
#' @param single_chain_list list passed from pmmh
#' @return list
#' @export

pmmh_single_chain <- function(single_chain_list) {

  ## Unlist object into local environment
  list2env(single_chain_list, .GlobalEnv)

  for(i in 2:ndraws){
    if(i == 2) {
      start_time <- Sys.time()
      param_values[1, ] <- ifelse(param_names %in% c("sigma", "phi"), 1, 0)
      ll_obj <- get_post_estim(data = data_list,
                               param_vec = param_values[1, ],
                               subj_index = 0,
                               param_name = "all",
                               nsim = nsim)

      pmll_old <- ll_obj$pmll_vec
      ll_old <- ll_obj$post_value
    }

    ## Propose new values
    old_vec <- param_values[i - 1, ]
    candidate_vec <- rnorm(num_param, old_vec, mh_prop_sd)
    accept_vec <- rep(0, num_param)
    for(j in 1:ncol(param_values)) {
      prop_vec <- old_vec
      prop_vec[j] <- candidate_vec[j]

      ## Get values for MH, omit of non-valid proposal
      if(param_names[j] %in% c("sigma", "phi") & prop_vec[j] <= 0) {
        ll_diff <- -1e12
      } else {
        #if(!param_names[j] %in% c("sigma", "phi", "alpha", "alpha_h")) prop_vec[j] <- 0
        ll_prop_obj <- get_post_estim(data = data_list,
                                      param_vec = prop_vec,
                                      param_name = param_names[j],
                                      subj_index = num_param_index[j],
                                      nsim = nsim,
                                      pmll_vec = pmll_old)

        #ll_old <- get_monte_estimate(sigma_old, phi_old)
        ll_diff <- ll_prop_obj$post_value - ll_old
      }

      ## MH acceptance step
      if(log(runif(1)) < ll_diff){
        old_vec[j] <- prop_vec[j]
        pmll_old <- ll_prop_obj$pmll_vec
        ll_old <- ll_prop_obj$post_value
        accept_vec[j] <- 1
      }
    }
    param_values[i, ] <- old_vec
    accept_values[i, ] <- accept_vec

    if(i %% 50 == 0) cat("Samples so far = ", i, "/", ndraws, "\n")
  }

  ## Create matrix with all draws
  if (keep_burn_in) {
    list(all_draws = param_values,
         all_accept = accept_values,
         eval_time = difftime(Sys.time(), start_time, units = "secs"))
  } else {
    list(all_draws = param_values[(burn_in + 1):ndraws, ],
         all_accept = accept_values[(burn_in + 1):ndraws, ],
         eval_time = difftime(Sys.time(), start_time, units = "secs"))
  }
}
