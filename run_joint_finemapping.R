#' Run Joint Fine-Mapping for All Modules
#'
#' @description
#' This function applies joint fine-mapping (`joint.finemapping`) to each module
#' defined in a `module.master.table`, and returns results for all modules.
#'
#' @param ann.data.obj A `SingleCellExperiment` or similar object containing credible set-level
#'        annotations and assays including `X`, `beta`, and `se`.
#' @param module.master.table A data.table or data.frame containing `cs_name` and `module` columns.
#'
#' @return A named list of results from `joint.finemapping()` for each module.
#'
#' @examples
#' \dontrun{
#' joint.results <- run_joint_finemapping_all_modules(
#'   ann.data.obj = ann.data,
#'   module.master.table = master_table
#' )
#'
#' # View most likely SNPs for module 1
#' joint.results[["1"]]$most.likely.snp.info
#' }
#'
#' @seealso \code{\link{joint.finemapping}}
#'
#' @export
#' 
run_joint_finemapping <- function(
    ann.data.obj,
    module.master.table
) {
  module.list <- split(module.master.table, by = "module", keep.by = FALSE)
  
  results <- lapply(module.list, function(module_table) {
    joint.finemapping(
      ann.data.obj = ann.data.obj,
      module.master.table = module_table
    )
  })
  
  names(results) <- names(module.list)
  
  ##### ADD HERE PART WHERE YOU ACTUALLY BIND MOST LIKELY SNPs and BETAS TO MASTER TABLE
  
  return(results)
}

### WOULD IT BE BETTER TO DO THE EXTRACTION FROM ANNDATA OUTSIDE OF THE LOOP? ALSO FOR RUN_COLOC_PERCOLATION
