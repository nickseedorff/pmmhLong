#' Calculate the gradient
#'
#' @param y vector, count for a specific user and time
#' @param f_cur vector, current estimate of the mode of the latent GP
#' @param inv_cov matrix, prior for GP precision matrix
#' @param mu vector, linear predictor without the GP
#'
#' @return vector
#' @export
#' @examples
#' inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
#' calc_gradient(c(1, 2), c(2, 3), inv_cov, c(-5, -4))

calc_gradient <- function(y, f_cur, inv_cov, mu) {
  as.vector(y - exp(f_cur + mu) - inv_cov %*% f_cur)
}

#' Caluate the hessian
#'
#' @param f_cur vector, current estimate of the mode of the latent GP
#' @param inv_cov matrix, prior for GP precision matrix
#' @param mu vector, linear predictor without the GP
#'
#' @return matrix
#' @export
#' @examples
#' inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
#' hessian(c(2, 3), inv_cov, c(-5, -4))

calc_hessian <- function(f_cur, inv_cov, mu) {
  -diag(exp(f_cur + mu)) - inv_cov
}

#' Take a single step in the newton raphson algorithm
#'
#' @param y vector, count for a specific user and time
#' @param f_cur vector, current estimate of the mode of the latent GP
#' @param inv_cov matrix, prior for GP precision matrix
#' @param mu vector, linear predictor without the GP
#' @param step_halve logical, step halving
#'
#' @return vector
#' @export
#' @examples
#' inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
#' single_step(c(1, 2), c(2, 3), inv_cov, c(-5, -4))

single_step <- function(y, f_cur, inv_cov, mu, step_halve = FALSE) {
  grad <- calc_gradient(y, f_cur, inv_cov, mu)
  inv_hess <- qr.solve(calc_hessian(f_cur, inv_cov, mu))

  if (step_halve) {
    as.vector(f_cur - inv_hess %*% grad / 2)
  } else {
    as.vector(f_cur - inv_hess %*% grad)
  }
}

#' Get posterior modes using laplace approximation
#'
#' @param y_mat matrix of count values where columns are the multivariate observations
#' @param cov_mat prior precision matrix for the latent GP
#' @param mean_mat matrix of linear predictors without the latent GP values
#'
#' @return vector
#' @export
#' @examples
#' y <- matrix(1:4, ncol = 2)
#' inv_cov <- qr.solve(matrix(c(1, 0.5, 0.5, 1), ncol = 2))
#' mean_mat <- y - 0.5
#' get_posterior_modes(y, inv_cov, mean_mat)

get_posterior_modes <- function(y_mat, inv_cov_mat, mean_mat) {

  ## Length of each GP component
  length_gp <- nrow(y_mat)
  num_obs <- ncol(y_mat)

  ## Store for posterior modes and
  post_f_modes <- matrix(0, nrow = nrow(y_mat), ncol = ncol(y_mat))
  prec_chol_array <- prec_array <- array(NA, dim = c(length_gp, length_gp,
                                                     ncol(y_mat)))

  ## Loop over all columns
  for(j in 1:num_obs){
    grad_norm <- 10
    i <- 0

    ## Observed values and linear preds
    y_vals <- y_mat[, j]
    mean_vals <- mean_mat[, j]
    fnew <- rep(0, length_gp)

    ## Iterate until converance
    while (i <= 50 & grad_norm > 1e-12){
      ## Store prior guess for comparison and update
      fold <- fnew
      fnew <- single_step(y_vals, fold, inv_cov_mat, mean_vals)

      lin_pred_old <- mean_vals + fold
      lin_pred_new <- mean_vals + fnew

      like_old <- sum(dpois(y_vals, exp(lin_pred_old), log = T)) - crossprod(lin_pred_old, inv_cov_mat) %*% lin_pred_old / 2
      like_new <- sum(dpois(y_vals, exp(lin_pred_new), log = T)) - crossprod(lin_pred_new, inv_cov_mat) %*% lin_pred_new / 2

      if (like_old > like_new) {
        fnew <- single_step(y_vals, fold, inv_cov_mat, mean_vals, step_halve = T)
      }

      ## Test convergence, update counter
      grad_norm <- crossprod(fold - fnew)
      i <- i + 1
    }

    ## Store posterior modes
    post_f_modes[, j] <- fnew
    prec_array[,, j] <- -calc_hessian(fnew, inv_cov_mat, mean_vals)
    prec_chol_array[,, j] <- chol(prec_array[,, j])
  }

  ## Return list of posterior modes and covariance matrices
  list(post_modes = post_f_modes, post_prec = prec_array,
       post_prec_chol = prec_chol_array)
}
