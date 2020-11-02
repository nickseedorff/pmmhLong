#' Calculate the gradient
#'
#' @param y vector, count for a specific user and time
#' @param f_guess vector, prior estimate of the mode of the latent GP
#' @param inv_cov matrix, current prior for GP precision matrix
#' @param mu vector, linear predictor without the GP
#'
#' @return vector
#' @export
#' @examples
#' inv_cov <- qr.solve(matrix(c(3, 1, 1, 3), ncol = 2))
#' calc_gradient(c(1, 2), c(2, 3), cov_m, c(-5, -4))

calc_gradient <- function(y, f_guess, inv_cov, mu) {
  as.vector(y - exp(mu + f_guess) - inv_cov %*% f_guess)
}
