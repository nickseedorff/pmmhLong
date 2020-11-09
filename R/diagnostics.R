#' Attain diagnostics for all chaings
#'
#' @param pmmh_list list returned from pmmh
#' @return list
#' @export

pmmh_diags <- function(pmmh_list) {
  pmmh_draws <- pmmh_list$all_draws
  pmmh_accepts <- pmmh_list$all_accept

  if (length(dim(pmmh_draws)) == 2) {
    diag_each_chain(pmmh_draws, pmmh_accepts,
                    pmmh_list$was_burn_in_kept,
                    pmmh_list$burn_in)
  } else {
    lapply(1:dim(pmmh_draws)[2], function(x) {
      draws_use <- pmmh_draws[, x, ]
      accepts_use <- pmmh_accepts[, x, ]
      diag_each_chain(draws_use, accepts_use,
                      pmmh_list$was_burn_in_kept,
                      pmmh_list$burn_in)
    })
  }
}

#' Attain diagnostics for each chain
#'
#' @param draws matrix of pmmh draws
#' @param accepts matrix of acceptances
#' @param was_burn_in_kept logical, were burn in values
#' @param burn_in burn in
#' @return list
#' @import coda
#' @export

diag_each_chain <- function(draws, accepts, was_burn_in_kept, burn_in) {

  ## Remove burn in if not previously done
  if (was_burn_in_kept) {
    draws <- draws[(burn_in + 1):nrow(draws), ]
    accepts <- accepts[(burn_in + 1):nrow(draws), ]
  }

  ## Get summary statistics
  summ_stat <- t(apply(draws, 2, function(x) {
    c(mean(x),
      as.numeric(quantile(x, probs = 0.025)),
      as.numeric(quantile(x, probs = 0.975)),
      sd(x))
  }))

  colnames(summ_stat) <- c("mean", "low_ci", "up_ci", "sd")

  ## Acceptance ratios
  accept_vec <- colMeans(accepts)

  ## convert to mcmc object, geweke diagnostic and correlations
  mc_obj <- mcmc(draws)
  neff <- effectiveSize(draws)
  gew_diag <- geweke.diag(mc_obj)$z
  auto_cor <- t(autocorr.diag(mc_obj, lag = c(1, 5)))
  cross_no_diag <- abs(crosscorr(mc_obj) - diag(rep(1, ncol(draws))))
  cross_cor_var <- apply(cross_no_diag, 2,
                         function(x) colnames(draws)[x == max(x)])
  cross_cor_val <- apply(cross_no_diag, 2, max)

  ## Add to dataframe and alter column names
  df <- data.frame(parameter = rownames(summ_stat),
                   summ_stat, neff = neff,
                   accept_prop = accept_vec, geweke_diag = gew_diag,
                   ac = auto_cor, cross_cor_var, cross_cor_val, row.names = NULL)
  colnames(df) <- gsub(pattern = "\\.", "_", tolower(colnames(df)))
  df
}
