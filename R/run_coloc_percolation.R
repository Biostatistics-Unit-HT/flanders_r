#' Run Coloc Percolation
#'
#' @description
#' This function groups traits by a `module` column and, for each group, extracts
#' the corresponding credible sets from annData/SCE object and runs coloc
#' percolation using a liberal `pp.H4.threshold` value.
#'
#' @param traits A data frame or tibble containing at least `module` and `cs_name` columns.
#' @param ad A `SingleCellExperiment` object containing annotation and assay data.
#' @param pp.H4.threshold Numeric threshold for colocalization (default = 0.5).
#' @param icoloc.reg Logical indicating whether to apply coloc regional filter (default = TRUE).
#'
#' @return A list of coloc percolation results, one per module.
#'
#' @examples
#' \dontrun{
#' results <- run_coloc_percolation(traits, ad)
#' }
#'
run_coloc_percolation <- function(traits,
                                  ad,
                                  pp.H4.threshold = 0.5,
                                  icoloc.reg = TRUE) {
  
  trait_groups <- dplyr::group_split(traits, traits$module)
  
  results <- lapply(trait_groups, function(x) {
    ad_sub <- extract.ann.data.labf.matrix(cs.list = x$cs_name, sce.obj = ad)
    coloc_percolation_matrix(ad_sub,
                             pp.H4.threshold = pp.H4.threshold,
                             icoloc.reg = icoloc.reg)
  })
  
  return(results)
}
