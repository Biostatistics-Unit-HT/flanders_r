#' Pad a Matrix with Missing Columns
#'
#' @description
#' Ensures that a matrix has all specified columns by adding any missing columns
#' and filling them with zeros (`0`). The function is useful for harmonizing
#' matrices prior to combining them (e.g., via `rbind`) when the column sets differ.
#'
#' @param mat A numeric matrix to be padded.
#' @param all_cols A character vector of all column names that should be present in the matrix.
#'
#' @return A matrix with the same number of rows as `mat` and exactly the column names
#'         specified in `all_cols`. Any columns originally missing from `mat` are added
#'         and filled with `0`. Columns are returned in the same order as `all_cols`.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("A", "B", "C")))
#' all_cols <- c("A", "B", "C", "D")
#' padded <- pad_matrix(mat, all_cols)
#' colnames(padded)  # Should be A, B, C, D
#' }
pad_matrix <- function(mat, all_cols) {
  missing_cols = setdiff(all_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    pad = matrix(0, nrow = nrow(mat), ncol = length(missing_cols),
                 dimnames = list(rownames(mat), missing_cols))
    mat = cbind(mat, pad)
  }
  mat = mat[, all_cols, drop = FALSE]  # Order columns consistently
  return(mat)
}

pad_matrix <- function(mat, all_cols) {
  missing_cols = setdiff(all_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    pad = matrix(0, nrow = nrow(mat), ncol = length(missing_cols),
                 dimnames = list(rownames(mat), missing_cols))
    mat = cbind(mat, pad)
  }
  mat = mat[, all_cols, drop = FALSE]  # Order columns consistently
  return(mat)
}