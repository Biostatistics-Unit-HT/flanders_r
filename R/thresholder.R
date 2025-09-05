#' Objective Function for Threshold Optimization
#'
#' Internal helper used by \code{\link{find_threshold}}. 
#' Given a numeric threshold and a vector of SNP scores, it splits the vector 
#' into two groups and returns the variance of residuals from a simple linear model.
#'
#' @param threshold Numeric. Threshold value used to split the vector into groups.
#' @param vector Numeric vector of SNP scores (e.g. lABFs).
#'
#' @return A single numeric value: the variance of residuals from the fitted model.
#'
#' @note This function is not exported. It is used internally by \code{find_threshold()}.
#'
#' @keywords internal
thresholder <- function(threshold, vector) {
  var(residuals(lm(vector ~ ifelse(vector > threshold, "A", "B"))))
}