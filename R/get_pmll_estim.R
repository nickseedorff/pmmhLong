#' Get estimates of the marginal likelihood values
#'
#' @param gen_data_mat_obj object with list of relevant matrices
#' @param nsim scalar, number of importance samples
#' @return vector
#' @export
#' @examples
#' #' get_pmll_estim(data_object, nsim = 20)

get_pmll_estim <- function(get_data_mat_obj, nsim = 500) {

  ## Get relevant matrices
  y_mat <- get_data_mat_obj$y_matrix
  cov_mat <- get_data_mat_obj$cov_matrix
  inv_cov_mat<- qr.solve(cov_mat)
  mean_mat <- get_data_mat_obj$mean_matrix
  log_fact_vec <- get_data_mat_obj$log_fact_vec

  ## Get posterior modes
  post_params <- get_posterior_modes(y_mat, inv_cov_mat, mean_mat)

  ## Qualifies for
  length_gp <- nrow(y_mat)
  num_obs <- ncol(y_mat)
  zero_vec <- rep(0, length_gp)

  ## Loop over each patient
  pmll_est <- sapply(1:num_obs, function(x){
    ## Posterior values of mode and precision matrix, determinant
    post_f <- post_params$post_modes[, x]
    post_prec <- post_params$post_prec[, , x]
    post_chol <- post_params$post_prec_chol[, , x]
    y_vals <- y_mat[, x]
    mean_vals <- mean_mat[, x]
    post_f_det <- prod(diag(post_chol))

    ## nsim simulations for each patient
    pmll_draws <- sapply(1:nsim, function(y){
      ## Get GP draws form linear prediction
      gp_draws <- post_f + backsolve(post_prec, rnorm(length_gp))
      #gp_draws <- mvrnorm(1, mu = post_f, Sigma = post_sig)
      linear_pred <- gp_draws + mean_vals

      ## Evaluate likelihood and densities
      pois_val <- sum(dpois(y_vals, exp(linear_pred), log = T))
      prior_val <- - crossprod(gp_draws, inv_cov_mat) %*% gp_draws / 2

      ## Evaluate posterior of GP
      gp_diff <- gp_draws - post_f
      post_val <- -crossprod(gp_diff, post_prec) %*% gp_diff / 2
      pois_val + prior_val - post_val
    })
    mean(exp(pmll_draws)) / post_f_det
  })

  ## Add determinant values to get person and time specific pmll estimates
  -sum(log(diag(chol(cov_mat)))) + log(pmll_est)
}
