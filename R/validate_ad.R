##' Validate an AnnData or SingleCellExperiment object for flanders
##'
##' Runs a set of structural checks to ensure the input object matches the
##' schema expected by downstream flanders functions. Specifically, it verifies
##' required columns in `var` (SNP axis) and `obs` (credible set axis), warns if
##' row/column names do not match their ID columns, prints brief per-chromosome
##' summaries, and executes a small end-to-end smoke test via
##' [anndata2coloc_input()] and [anndata2coloc()].
##'
##' Accepted inputs are either a Python-backed AnnData object (via reticulate)
##' or a SingleCellExperiment created with zellkonverter. In both cases, the
##' conventions are:
##' - `obs` (rows) represent credible sets; must include: `cs_name`, `chr`,
##'   `start`, `end`, `study_id`, `phenotype_id`, `min_res_labf`, `coverage`,
##'   `logsum.logABF`.
##' - `var` (columns) represent SNPs; must include: `snp`, `chr`, `pos`.
##' - `rownames(obs) == obs$cs_name` and `rownames(var) == var$snp` are advised
##'   (the function warns if they differ).
##'
##' The function stops with an error if any required columns are missing or if
##' the smoke tests return empty outputs. Otherwise, it prints a readable
##' summary and invisibly returns a compact list of diagnostics.
##'
##' @param ad An AnnData or SingleCellExperiment object. For AnnData, required
##'   fields must be present in `ad$obs` and `ad$var`; for SingleCellExperiment,
##'   they must be present in `colData(ad)` and `rowData(ad)` respectively.
##'
##' @return Invisibly returns a list with elements:
##'   - `n_cs`, `n_snp`: counts of credible sets and SNPs
##'   - `cs_per_chr`, `snp_per_chr`: tables per chromosome
##'   - `coloc_input_rows`: number of rows returned by [anndata2coloc_input()]
##'   - `coloc_test_rows`: number of rows returned by [anndata2coloc()] on a
##'     10-row subset
##'   The function has side effects: it prints validation output and will `stop()`
##'   if checks fail.
##'
##' @seealso [anndata2coloc_input()], [anndata2coloc()],
##'   zellkonverter::readH5AD, anndata::read_h5ad
##'
##' @examples
##' \dontrun{
##' # Example 1: SingleCellExperiment from .h5ad
##' library(zellkonverter)
##' library(SingleCellExperiment)
##' sce <- readH5AD("/path/to/my_anndata.h5ad", reader = "R")
##' validate_anndata(sce)
##'
##' # Example 2: reticulate-backed AnnData
##' library(anndata)
##' ad <- read_h5ad("/path/to/my_anndata.h5ad")
##' validate_anndata(ad)
##' }
##'
##' @export
validate_anndata <- function(ad) {
  cat("---- Validating AnnData / SingleCellExperiment format ----\n")
  
  # detect type
  if (inherits(ad, "SingleCellExperiment")) {
    var <- as.data.frame(rowData(ad))
    obs <- as.data.frame(colData(ad))
    var_rownames <- rownames(ad)
    obs_rownames <- colnames(ad)
  } else if (inherits(ad, "AnnData")) {
    var <- as.data.frame(ad$var)
    obs <- as.data.frame(ad$obs)
    var_rownames <- rownames(ad$var)
    obs_rownames <- rownames(ad$obs)
  } else {
    stop("Unsupported object type: must be AnnData or SingleCellExperiment")
  }
  
  # ---- check required var ----
  required_var <- c("snp", "chr", "pos")
  missing_var <- setdiff(required_var, colnames(var))
  if (length(missing_var) > 0) {
    cat("ERROR: Missing columns in var:", paste(missing_var, collapse=", "), "\n")
    cat("Available columns:\n")
    print(colnames(var))
    cat("\nHead of var:\n")
    print(head(var, 5))
    stop("Missing required columns in var: ", paste(missing_var, collapse=", "))
  }
  
  # rownames vs snp (warn only)
  if (!all(var_rownames == var$snp)) {
    cat("WARNING: rownames(var) do not match var$snp\n")
    print(head(data.frame(rownames = var_rownames, snp = var$snp), 5))
  } else {
    cat("OK: var rownames match var$snp\n")
  }
  
  # ---- check required obs ----
  required_obs <- c("cs_name", "chr", "start", "end", 
                    "study_id", "phenotype_id", "min_res_labf", 
                    "coverage", "logsum.logABF")
  missing_obs <- setdiff(required_obs, colnames(obs))
  if (length(missing_obs) > 0) {
    cat("ERROR: Missing columns in obs:", paste(missing_obs, collapse=", "), "\n")
    cat("Available columns:\n")
    print(colnames(obs))
    cat("\nHead of obs:\n")
    print(head(obs, 5))
    stop("Missing required columns in obs: ", paste(missing_obs, collapse=", "))
  }
  
  # colnames vs cs_name (warn only)
  if (!all(obs_rownames == obs$cs_name)) {
    cat("WARNING: obs names do not match obs$cs_name\n")
    print(head(data.frame(obs_names = obs_rownames, cs_name = obs$cs_name), 5))
  } else {
    cat("OK: obs names match obs$cs_name\n")
  }
  
  # ---- summary ----
  cat("\n---- Summary ----\n")
  n_cs   <- nrow(obs)
  n_snp  <- nrow(var)
  
  cat("Credible sets (CS):", n_cs, "\n")
  cat("SNPs:", n_snp, "\n\n")
  
  cs_per_chr <- table(obs$chr)
  snp_per_chr <- table(var$chr)
  
  cat("CS per chromosome:\n")
  print(cs_per_chr)
  cat("\nSNPs per chromosome:\n")
  print(snp_per_chr)
  
  # ---- integration tests ----
  cat("\n---- Integration tests ----\n")
  coloc_input <- anndata2coloc_input(ad)
  if (!is.data.frame(coloc_input) || nrow(coloc_input) == 0) {
    stop("Integration test failed: anndata2coloc_input returned empty")
  } else {
    cat("anndata2coloc_input: OK (", nrow(coloc_input), "rows)\n")
  }
  
  test_input <- head(coloc_input, 10)
  coloc_res <- anndata2coloc(ad, test_input)
  if (!is.data.frame(coloc_res) || nrow(coloc_res) == 0) {
    stop("Integration test failed: anndata2coloc returned empty")
  } else {
    cat("anndata2coloc: OK (", nrow(coloc_res), "rows)\n")
  }
  
  cat("\nDONE: Validation completed\n")
  
  invisible(list(
    n_cs = n_cs,
    n_snp = n_snp,
    cs_per_chr = cs_per_chr,
    snp_per_chr = snp_per_chr,
    coloc_input_rows = nrow(coloc_input),
    coloc_test_rows = nrow(coloc_res)
  ))
}