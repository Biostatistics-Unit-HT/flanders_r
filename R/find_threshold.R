#' Identify Expanded Credible Set SNPs via Thresholding
#'
#' Given a vector of SNP scores (e.g. lABFs), this function determines an 
#' optimal threshold that best separates the values into two groups. 
#' SNPs above the threshold are considered part of the expanded credible set.
#'
#' @param vector Numeric vector of SNP-level scores.
#' @param original_idx Optional integer vector of SNP indices corresponding 
#'   to the original credible set. If no expansion is possible 
#'   (e.g. variance of \code{vector} is zero), these indices are returned instead.
#'
#' @return An integer vector of SNP indices corresponding to the expanded credible set.  
#'   If no expansion is possible, the original indices are returned.
#'
#' @seealso \code{\link{expand_cs}}
#' @export
#'
#' @examples
#' lbf_vec <- rnorm(100)
#' cs_idx <- which(lbf_vec > 2)
#' find_threshold(lbf_vec, original_idx = cs_idx)
find_threshold <- function(vector, original_idx = NULL) {
  if (var(vector) != 0) {
    thres <- optim(
      fn = function(th) thresholder(th, vector),
      p = 0,
      method = "Brent",
      lower = min(vector),
      upper = max(vector)
    )
    idx <- which(vector > thres$par)
    return(idx)
  } else {
    # Return original indices if no expansion is possible
    return(original_idx)
  }
}