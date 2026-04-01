#' Add PIP and is.cs99 layers to a SingleCellExperiment (read by zellkonverter::readH5AD)
#'
#' This function expects `sce` to be the object returned by `readH5AD(..., reader = "R")` and
#' that `colData(sce)$logsum.logABF` exists. It computes per-column PIP (by subtracting the
#' provided `logsum.logABF` from the column lABF values), stores the PIP matrix as
#' `sce@assays@data$PIP` and a logical 99% cumulative-PIP mask as `sce@assays@data$is.cs99`.
#'
#' @param sce A SingleCellExperiment object (from `readH5AD(..., reader = "R")`).
#' @param threshold Numeric coverage threshold for the credible set (default 0.99).
#' @return The input `sce` with two added assays/layers: `PIP` (numeric, same sparse format)
#'         and `is.cs99` (sparse logical matrix).
#' @useDynLib flanders, .registration = TRUE
#' @examples
#' \dontrun{
#' library(zellkonverter)
#' sce <- readH5AD("/path/to/my_anndata.h5ad", reader = "R")
#' sce <- anndata_add_pip_iscs99(sce, threshold = 0.99)
#' }
#' @export

anndata_add_pip_iscs99 <- function(sce, threshold = 0.99) {
  if (is.character(sce)) stop("`sce` must be a SingleCellExperiment object, not a path")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required")

  # Extract lABF matrix and per-column logsum.logABF vector
  X <- sce@assays@data$X
  if (is.null(X)) stop("sce@assays@data$X not found")

  v <- colData(sce)$logsum.logABF
  if (is.null(v)) stop("colData(sce)$logsum.logABF not found")
  if (length(v) != ncol(X)) warning("Length of logsum.logABF does not match ncol(X)")

  # Work on a copy of X to build PIP (log-space)
  PIP_log <- X
  PIP_log@x <- PIP_log@x - rep(as.numeric(v), diff(PIP_log@p))

  # Create logical sparse matrix for is.cs99
  nr <- nrow(PIP_log)
  nc <- ncol(PIP_log)
  out <- Matrix::Matrix(FALSE, nrow = nr, ncol = nc, sparse = TRUE)

  # Call optimized C++ implementation to compute kept row indices per column
  p <- as.integer(PIP_log@p)
  ii <- as.integer(PIP_log@i)
  xx <- as.numeric(PIP_log@x)

  res <- compute_cs99_cpp(p, ii, xx, as.integer(nc), as.numeric(threshold))

  # assemble sparse logical matrix from results
  if (length(res) == 0) {
    out <- Matrix::Matrix(FALSE, nrow = nr, ncol = nc, sparse = TRUE)
  } else {
    j_counts <- vapply(res, length, integer(1))
    if (sum(j_counts) == 0) {
      out <- Matrix::Matrix(FALSE, nrow = nr, ncol = nc, sparse = TRUE)
    } else {
      i_vals <- unlist(res, use.names = FALSE)
      j_vals <- rep(seq_along(res), times = j_counts)
      out <- Matrix::sparseMatrix(i = i_vals, j = j_vals, x = TRUE,
                                 dims = c(nr, nc), index1 = TRUE)
    }
  }

  # Assign layers to sce
  sce@assays@data$PIP <- PIP_log
  sce@assays@data$is.cs99 <- out

  # Convert PIP from log to linear scale in-place
  sce@assays@data$PIP@x <- exp(sce@assays@data$PIP@x)

  return(sce)
}

# Deprecated: previously a standalone script; use `anndata_add_pip_iscs99(sce, threshold)` instead.