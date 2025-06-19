#' Assign Unique Module IDs to Rows with Missing Module Information
#'
#' @description
#' This function replaces `NA` values in the `module` column of a data frame
#' with unique numeric IDs, starting from a specified value. It is useful for
#' assigning module labels to traits or regions that were not part of any
#' colocalized module (i.e., singletons).
#'
#' @param master.table_w_modules A `data.frame` or `data.table` with a column named `module`,
#'   where some entries may be `NA`.
#' @param start_value A numeric value indicating the starting number to assign to
#'   `NA` entries in the `module` column. The values assigned will increment by 1
#'   for each subsequent `NA`. If `length(na_indices) = 0`, no changes are made.
#'
#' @return A modified version of `master.table_w_modules` where all `NA` values
#'   in the `module` column are replaced with unique numeric module IDs.
#'
#' @examples
#' df <- data.frame(cs_name = c("a", "b", "c"), module = c(1, NA, NA))
#' replace_nas_with_modules(df, start_value = 100)
#'
#' @keywords internal
#'
replace_nas_with_modules <- function(master.table_w_modules, start_value=NULL) {
  na_indices <- which(is.na(master.table_w_modules[["module"]]))
  master.table_w_modules[["module"]][na_indices] <- seq(from = start_value, length.out = length(na_indices))
  return(master.table_w_modules)
}