#' Merge and filter meta-analysis CS into an existing SCE
#'
#' This function:
#' \enumerate{
#'   \item Identifies any CS that appear in more than one meta-analysis pairing, logs them, and excludes them.
#'   \item Drops from the original SCE all CS involved in the retained meta-analysis runs.
#'   \item Builds new assay matrices for labf (primary assay \code{X}), \code{beta}, and \code{se}.
#'   \item Constructs a new SCE for the meta-CS and column-binds it to the filtered SCE.
#' }
#'
#' @param results A list of `cs_metanalysis` outputs. Each element must be a list with components:
#'   \describe{
#'     \item{list_of_cs_metaanalyzed}{Character vector of two original CS names.}
#'     \item{cs.name}{The new meta-CS name (length 1).}
#'     \item{meta.results}{A data.frame with columns `SNPID`, `beta`, `se`, `logp`, `labf`, and `PPI`.}
#'   }
#' @param sce A `SingleCellExperiment` with:
#'   \itemize{
#'     \item `rownames(sce)` equal to SNPIDs
#'     \item `colnames(sce)` equal to original CS names
#'     \item assays named `"X"`, `"beta"`, and `"se"` present
#'   }
#' @return A `SingleCellExperiment` where:
#' \describe{
#'   \item{assays}{Original assays (`X`, `beta`, `se`) minus the dropped CS columns, plus new columns for each meta-CS.}
#'   \item{colData}{Extended to include one row per meta-CS, with `cs_name` and related metadata.}
#' }
#' @import SingleCellExperiment
#' @export

ma2anndata <- function(results, sce) {
  # 1) blacklist CS appearing >1× across all pairings
  all_cs <- unlist(lapply(results, function(x) x$list_of_cs_metaanalyzed))
  counts <- table(all_cs)
  bad_cs <- names(counts)[counts != 1]
  if (length(bad_cs) > 0) {
    message("Excluding CS that appear multiple times:",
            paste0("
 - ", bad_cs, collapse = ""))
    results <- Filter(function(x) !any(x$list_of_cs_metaanalyzed %in% bad_cs), results)
  }

  # 2) drop original CS columns from sce
  to_drop <- unique(unlist(lapply(results, `[[`, "list_of_cs_metaanalyzed")))
  sce_filtered <- sce[, !(colnames(sce) %in% to_drop)]

  # 3) prepare new meta-CS matrices
  snps    <- rownames(sce_filtered)
  meta_cs <- vapply(results, `[[`, "cs.name", FUN.VALUE = character(1))
  n_meta  <- length(meta_cs)
  labf_mat <- matrix(NA_real_, nrow = length(snps), ncol = n_meta,
                     dimnames = list(snps, meta_cs))
  beta_mat <- labf_mat
  se_mat   <- labf_mat

  for (i in seq_along(results)) {
    df  <- results[[i]]$meta.results
    sel <- match(df$SNPID, snps, nomatch = 0)
    labf_mat[sel, i] <- df$labf
    beta_mat[sel, i] <- df$beta
    se_mat[sel, i]   <- df$se
  }

  # 4) build a new SCE for meta-CS
  # replicate filtered colData structure
  orig_cd <- colData(sce_filtered)
  cd_names <- colnames(orig_cd)
  cd_list <- lapply(cd_names, function(x) rep(NA, n_meta))
  names(cd_list) <- cd_names
  sce2add_cd <- S4Vectors::DataFrame(cd_list, row.names = meta_cs)

  # fill known colData fields
  sce2add_cd$cs_name <- meta_cs
  for (i in seq_along(results)) {
    name_i <- meta_cs[i]
    res    <- results[[i]]
    sce2add_cd[name_i, "chr"]           <- as.character(res$chr)
    sce2add_cd[name_i, "start"]         <- res$start
    sce2add_cd[name_i, "end"]           <- res$end
    parts                                <- strsplit(name_i, "::", fixed = TRUE)[[1]]
    sce2add_cd[name_i, "study_id"]      <- parts[1]
    sce2add_cd[name_i, "phenotype_id"]  <- parts[2]
    sce2add_cd[name_i, "min_res_labf"]  <- res$min.abf
  }

  sce2add <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(X = labf_mat, beta = beta_mat, se = se_mat),
    rowData = rowData(sce_filtered),
    colData = sce2add_cd
  )

  # 5) combine filtered sce and new meta-CS SCE
  sce_out <- SingleCellExperiment::cbind(sce_filtered, sce2add)
  return(sce_out)
}
