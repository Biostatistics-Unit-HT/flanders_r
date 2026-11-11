#' Extract conditional effect estimates from a SuSiE fit
#'
#' Extracts posterior means and posterior standard deviations conditional on
#' inclusion for all variants in a selected SuSiE single-effect component.
#' Estimates are rescaled using the column scale factors stored in the fit.
#'
#' @param sus A SuSiE fitted model returned by \code{susieR::susie_rss()},
#'   containing \code{mu}, \code{mu2}, and \code{X_column_scale_factors}.
#' @param L_index A single integer identifying the single-effect component
#'   to extract, between 1 and \code{nrow(sus$mu)}. This is a row index in
#'   \code{sus$mu}, rather than a position in the list of credible sets.
#'
#' @details
#' For each variant in the selected component, the posterior mean is
#' calculated as \code{mu / X_column_scale_factors}, and the posterior
#' standard deviation as
#' \code{sqrt(mu2 - mu^2) / X_column_scale_factors}.
#'
#' Both quantities are conditional on the variant being included in the
#' selected component; they are not weighted by posterior inclusion
#' probabilities (\code{sus$alpha}).
#'
#' The output scale depends on the inputs supplied to
#' \code{susieR::susie_rss()}; rescaling does not necessarily recover
#' effect sizes on the original trait and genotype scales.
#'
#' @return A named list containing two numeric vectors, each with one entry
#'   per variant in the same order as the columns of \code{sus$mu}:
#'   \describe{
#'     \item{beta}{Posterior means conditional on inclusion, after rescaling.}
#'     \item{se}{Posterior standard deviations conditional on inclusion,
#'       after rescaling.}
#'   }
#'
#' @seealso \code{\link[susieR]{susie_rss}}
#'
#' @export
get_beta_se_susie <- function(sus, L_index){
  
  se = sqrt(sus$mu2[L_index,] - (sus$mu[L_index,])^2)/sus$X_column_scale_factors
  beta = (sus$mu[L_index,]/sus$X_column_scale_factors)
  
  return(list(beta = beta, se = se))
  
}