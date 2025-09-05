#' Objective Function for Threshold Optimization
#'
#' Internal helper function used in credible set expansion.  
#' Given a threshold, it fits a simple linear model splitting the input vector 
#' into two groups and returns the variance of the residuals. This is used as 
#' the objective function for \code{optim()} in \code{find_threshold()}.
#'
#' @param threshold Numeric. Threshold value used to split the vector into groups.
#'
#' @return A numeric value: the variance of residuals from the fitted model.
#'
#' @note This function is not meant to be called directly by the user.  
#' It is called internally by \code{\link{find_threshold}}.
#'
#' @keywords internal
thresholder <- function(threshold) {
    var(residuals(lm(vector ~ ifelse(vector > threshold, "A", "B"))))
}