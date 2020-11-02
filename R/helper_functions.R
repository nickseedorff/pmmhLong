#' Build data matrices for ingestion by gen_pmll_estim()
#'
#' @param data_obj list of objects to make relevant data matrices
#' @param sigma scalar
#' @param phi scalar
#' @param mean_vec vector, linear predictor without the GP
#'
#' @return list of matrices
#' @export

gen_data_mat <- function(data_obj, sigma, phi, mean_vec) {
  length_gp <- nrow(data_obj$dist_mat)
  y_matrix <- matrix(data_obj$Y, nrow = length_gp)
  cov_matrix <- sigma * exp(-data_obj$dist_mat ^ 2/ phi)
  mean_matrix <- matrix(mean_vec, nrow = length_gp, ncol = ncol(y_matrix))
  log_fact_vec <- apply(y_matrix, 2, function(x) sum(lfactorial(x)))
  list(y_matrix = y_matrix, cov_matrix = cov_matrix, mean_matrix = mean_matrix,
       log_fact_vec = log_fact_vec)
}
