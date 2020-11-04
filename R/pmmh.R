#' Get values for MH accept step
#'
#' @param data data.frame or matrix with described column names
#' @param distance_mat matrix, distance between multivariate outcomes
#' @param outcome char, name of the outcome variable
#' @param subject_index char, name of the column denoting subjects
#' @param time_index char, name of the column denoting time
#' @param position_index char, name of the column denoting position
#' @param ndraws number of draws
#' @param nsim number of importance samples
#' @param offset_term char, name of column that is used as an offset
#' @return list
#' @export

pmmh <- function(data, distance_mat,
                 outcome = "y", subject_index = "patient",
                 time_index = "time", position_index = "position",
                 ndraws = 200, nsim = 10, offset_term = NULL) {

  ## Scalars needed later and prepare list of relevant objects
  length_gp <- nrow(distance_mat)
  num_subjects <- length(unique(data[, subject_index]))
  data_list <- prepare_data(data, distance_mat, length_gp, num_subjects,
                            outcome, subject_index, time_index,
                            position_index, offset_term)

  ## Store draws from posterior
  param_names <- data_list$param_names
  num_param <- length(param_names)
  num_param_index <- rep(0, num_param)
  num_param_index[1:(2 * num_subjects)] <- rep(1:num_subjects, 2)
  accept_values <- param_values <- matrix(NA, nrow = ndraws, ncol = num_param)

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
    candidate_vec <- rnorm(num_param, old_vec,
                           c(rep(0.5, 2 * num_subjects), 0.1, 0.5, 2, 0.5))
    accept_vec <- rep(0, num_param)
    for(j in 1:ncol(param_values)) {
      prop_vec <- old_vec
      prop_vec[j] <- candidate_vec[j]

      ## Get values for MH, omit of non-valid proposal
      if(param_names[j] %in% c("sigma", "phi") & prop_vec[j] <= 0) {
        ll_diff <- -100000000
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

    if(i %% 25 == 0) print(i)
  }

  ## Create matrix with all draws
  param_names <- c(paste0("alpha", 1:num_subjects),
                   paste0("beta", 1:num_subjects),
                   "sigma", "phi", "alpha_h")

  list(param_names = param_names,
       all_draws = param_values,
       all_accept = accept_values,
       eval_time = difftime(Sys.time(), start_time, units = "mins"))
}
