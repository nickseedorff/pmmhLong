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
#' @import dplyr magrittr tidyr ggplot2
#' @export


pmmh_trace <- function(pmmh_list, parameter = "alpha") {
  dat <- gen_data(sigma = 0.8, phi = 2, seed = 20,
                  patients = 10, time_points_per_patient = 4,
                  alpha_hier = 0.5, beta_hier = 0.3)

  pmmh_list <- pmmh(dat$all_dat, dat$distance_matrix, chains = 3, ndraws = 6000, burn_in = 3000,
                    nsim = 25)

  pmmh_diags(pmmh_list)

  color_scheme_set("brewer-Spectral")
  mcmc_intervals(pmmh_list$all_draws)
  mcmc_trace(pmmh_list$all_draws, regex_pars = "alpha.")
  mcmc_trace(pmmh_list$all_draws, regex_pars = "beta.")
  mcmc_acf(pmmh_list$all_draws, regex_pars = "alpha.")
  mcmc_acf(pmmh_list$all_draws, regex_pars = "beta.")
  mcmc_trace(pmmh_list$all_draws, regex_pars = "sigma")
  mcmc_trace(pmmh_list$all_draws, regex_pars = "phi")


  dat <- gen_data(sigma = 0.6, phi = 2, seed = 4,
                  patients = 10, time_points_per_patient = 3,
                  alpha_hier = 0.3, beta_hier = 0.2)

  pmmh_list <- pmmh(dat$all_dat, dat$distance_matrix, nsim = 25, ndraws = 5000)
  pmmh_list2 <- pmmh(dat$all_dat, dat$distance_matrix, nsim = 25, ndraws = 5000)

  ## convert to dataframe
  draws <- pmmh_list$all_draws
  colnames(draws) <- pmmh_list$param_names
  draws <- draws[, colnames(draws) == parameter, drop = FALSE]

  ## Update names
  if (ncol(draws) > 1) {
    colnames(draws) <- paste0(parameter, 1:ncol(draws))
  }

  draws2 <- pmmh_list2$all_draws
  colnames(draws2) <- pmmh_list$param_names
  draws2 <- draws2[, colnames(draws2) == parameter, drop = FALSE]

  ## Update names
  if (ncol(draws2) > 1) {
    colnames(draws2) <- paste0(parameter, 1:ncol(draws2))
  }

  arr <- abind::abind(list(draws, draws2))

  arr <- array(NA, dim = c(nrow(draws), 2, ncol(draws)), dimnames = list(iterations = NULL,
                                                           chains = c("one", "two"),
                                                           parameters = colnames(draws)))

  dimnames(arr)
  dim(arr)
  arr[,1,] <- draws
  arr[,2,] <- draws2
  dim(arr)

  ## Bayes plot can use an array
  mcmc_areas(draws)
  mcmc_trace(draws)
  mcmc_intervals(draws)
  mcmc_hist(draws)

  color_scheme_set("brewer-Spectral")
  mcmc_trace(arr, pars = c("sigma", "phi", "alpha_h", "beta_h"))

  mcmc_areas(arr)
  mcmc_trace(arr, pars = c("sigma", "phi", "alpha_h", "beta_h"))
  mcmc_intervals(draws)
  mcmc_hist(draws)

  mc <- mcmc(draws[501:5000, ])
  geweke.diag(mc)

  mc2 <- mcmc(draws2[501:5000, ])
  mc_list <- as.mcmc.list(list(mc, mc2))
  gelman.diag(mc_list)


  ## Add columns and pivot
  draws$iteration <- 1:nrow(draws)
  plot_df <- draws %>%
    pivot_longer(cols = starts_with(parameter))

  ## Plot
  ggplot(plot_df, aes(x = iteration, y = value)) +
    geom_line() +
    facet_grid(name ~ ., scales = "free_y")
}

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
#' @import dplyr magrittr tidyr ggplot2
#' @export


pmmh_hist <- function(pmmh_list, parameter = "alpha") {

  ## convert to dataframe
  draws <- as.data.frame(pmmh_list$all_draws)
  colnames(draws) <- pmmh_list$param_names
  draws <- draws[, colnames(draws) == parameter, drop = FALSE]

  ## Update names
  if (ncol(draws) > 1) {
    colnames(draws) <- paste0(parameter, 1:ncol(draws))
  }

  ## Add columns and pivot
  draws$iteration <- 1:nrow(draws)
  plot_df <- draws %>%
    pivot_longer(cols = starts_with(parameter))

  ## Plot
  ggplot(plot_df, aes(value)) +
    geom_histogram(bins = 100) +
    facet_grid(name ~ ., scales = "free_y")
}

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
#' @import dplyr magrittr tidyr ggplot2
#' @export


pmmh_geweke <- function(pmmh_list, parameter = "alpha", start = 101) {

  ## convert to mcmc object for coda
  draws <- pmmh_list$all_draws
  colnames(draws) <- pmmh_list$param_names
  mc <- mcmc(draws, start = 250)
  geweke.diag(mc)


  mcmc_areas(draws)

  autocorr.plot(mc)

  ## Need to create multiple chains for additional diagnostic plots

  traceplot(mc)
  draws <- draws[, colnames(draws) == parameter, drop = FALSE]

  ## Update names
  if (ncol(draws) > 1) {
    colnames(draws) <- paste0(parameter, 1:ncol(draws))
  }

  ## Add columns and pivot
  draws$iteration <- 1:nrow(draws)
  plot_df <- draws %>%
    top_frac()
    pivot_longer(cols = starts_with(parameter))

  ## Plot
  ggplot(plot_df, aes(value)) +
    geom_histogram(bins = 100) +
    facet_grid(name ~ ., scales = "free_y")
}


library(bayesplot)
