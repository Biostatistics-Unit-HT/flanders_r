#' Run Coloc Percolation
#'
#' @description
#' Groups traits by a `module` column and runs colocalization percolation for each group.
#' Internally, the function first extracts log-approximated Bayes factors (lABFs) for all
#' credible sets across all modules from the input `SingleCellExperiment` object (`ad`),
#' and then subsets the matrix per group. This avoids redundant data extraction in loops.
#'
#' For each module, the corresponding lABF subset is passed to `coloc_percolation_matrix()`,
#' which constructs pairwise coloc matrices, applies optional regional regularization, and
#' identifies colocalization cliques.
#'
#' @param traits A `data.frame` or `tibble` containing at least the columns `module` and `cs_name`.
#'        Each row defines a credible set and its assigned module.
#' @param ad A `SingleCellExperiment` object containing credible set-level annotations,
#'        and assays including `X` (lABF), `beta`, and `se`.
#' @param pp.H4.threshold Numeric threshold for considering colocalization via posterior probability H4.
#'        Default is 0.5.
#' @param icoloc.reg Logical flag indicating whether to apply coloc regional regularization to each
#'        ABF vector. Default is TRUE.
#'
#' @return A list of coloc percolation result objects, one per module. Each includes a colocalization
#'         matrix, detected cliques, assigned trait groups, and ABF submatrices.
#'
#' @seealso \code{\link{coloc_percolation_matrix}}, \code{\link{extract.ann.data.labf.matrix}}
#'
#' @examples
#' \dontrun{
#' results <- run_coloc_percolation(traits = trait_table, ad = ann_data)
#' results[[1]]$coloc.matrix
#' }
#'
#' @export
#'
run_coloc_percolation <- function(traits,
                                  ad,
                                  pp.H4.threshold = 0.5,
                                  icoloc.reg = TRUE) {
  
  trait_groups <- dplyr::group_split(traits, traits$module)
  
  # Extract lABFs for all cs_names at once
  all_cs <- traits$cs_name
  full_labf_matrix <- extract.ann.data.labf.matrix(cs.list = all_cs, sce.obj = ad)
  
  # Loop coloc percolation by modules on the precomputed matrix
  results <- lapply(trait_groups, function(x) {
    
    cs_names <- x$cs_name
    subset_labf_matrix <- full_labf_matrix[, cs_names, drop = FALSE]
    
    coloc_percolation_matrix(
      finemap.data = subset_labf_matrix,
      pp.H4.threshold = pp.H4.threshold,
      icoloc.reg = icoloc.reg)
  })
  
  return(results)
}
