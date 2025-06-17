#' Convert Log-Bayes Factors to Posterior Probabilities
#'
#' @description
#' This function transforms a vector of log-approximated Bayes factors (lABFs)
#' into a vector of posterior probabilities. It normalizes the lABFs by
#' dividing each `exp(lABF)` by the sum of `exp(all lABFs)`, performing these
#' operations in log-space for numerical stability.
#'
#' @param bf A numeric vector of log-approximated Bayes factors (lABFs).
#'   These lABFs represent the evidence for different hypotheses
#'   (e.g., colocalization hypotheses).
#'
#' @return A numeric vector of posterior probabilities, corresponding to the
#'   input `bf` values. The sum of all elements in the returned vector will
#'   be approximately 1.
#'
#' @details
#' To convert Bayes factors (BFs) to posterior probabilities (PPs), assuming
#' equal prior probabilities for all hypotheses, the formula is:
#' `PP_i = BF_i / sum(BF_j for all j)`
#'
#' When working with log-Bayes factors (lABFs), this translates to:
#' `log(PP_i) = log(BF_i) - log(sum(BF_j))`
#' `log(PP_i) = log(exp(lABF_i)) - log(sum(exp(lABF_j)))`
#' `log(PP_i) = lABF_i - log(sum(exp(lABF_j)))`
#'
#' This function implements the numerically stable calculation of `log(sum(exp(lABF_j)))`
#' using the `coloc::logsum` function, and then computes `exp(lABF_i - log(sum(exp(lABF_j))))`
#' to obtain the posterior probabilities. This approach prevents issues like
#' underflow or overflow that can occur when directly exponentiating and summing
#' very small or very large log-values.
#'
#' @seealso
#' \code{\link[coloc]{logsum}} (for the underlying numerically stable log-sum-exp calculation)
#'
#' @export
#'
logbf_to_pp <- function(bf=all.coloc.join$ABF.tot) {
  denom=coloc:::logsum(bf)
  exp(bf  - denom)
}

