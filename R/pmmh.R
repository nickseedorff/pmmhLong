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
                 ndraws = 200, burn_in = 51, nsim = 10, offset_term = NULL,
                 chains = 1, keep_burn_in = FALSE) {

  if (burn_in >= ndraws) {
    stop("burn_in must be smaller than ndraws")
  }

  ## Scalars needed later and prepare list of relevant objects
  length_gp <- nrow(distance_mat)
  num_subjects <- length(unique(data[, subject_index]))
  data_list <- prepare_data(data, distance_mat, length_gp, num_subjects,
                            outcome, subject_index, time_index,
                            position_index, offset_term)
  single_chain_list <- prepare_storage(data_list, ndraws, num_subjects, nsim,
                                       burn_in, keep_burn_in)

  if (chains == 1) {
    results_list <- pmmh_single_chain(single_chain_list)

  } else if (chains > 1){
    tmp_list <- lapply(1:chains, function(x) {
      pmmh_single_chain(single_chain_list)
    })

    ## Convert to arrays with appropriate names
    draws1 <- tmp_list[[1]]$all_draws
    all_draws <- array(NA, dim = c(nrow(draws1), chains, ncol(draws1)),
                       dimnames = list(iterations = NULL,
                                       chains = paste0("chain", 1:chains),
                                       paramters = colnames(draws1)))
    all_accept <- all_draws
    time_vec <- vector(length = chains)
    for(i in 1:chains) {
      all_draws[, i,] <- tmp_list[[i]]$all_draws
      all_accept[, i,] <- tmp_list[[i]]$all_accept
      time_vec[i] <- paste0("Chain ", i, " runtime = ",
                            round(as.numeric(tmp_list[[i]]$eval_time), 1),
                            " seconds")
    }
    results_list = list(all_draws = all_draws, all_accept = all_accept,
                        runtime = time_vec)
  }
  ## Add burn_in info
  burn_info <- list(burn_in = single_chain_list$burn_in,
                    was_burn_in_kept = single_chain_list$keep_burn_in)

  append(results_list, burn_info)
}
