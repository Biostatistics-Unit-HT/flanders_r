#' Run Joint Fine-Mapping and Aggregate Results
#'
#' @description
#' Applies joint fine-mapping (`joint.finemapping`) to each module defined in `module.master.table`.
#' The function merges most likely SNP information back into the input table, and constructs a sparse matrix
#' of joint log10 Bayes Factors (`lABF`) for all credible sets across modules.
#'
#' @param ann.data.obj A `SingleCellExperiment` or similar object containing credible set-level
#'        annotations and assays, including `X`, `beta`, and `se`.
#' @param module.master.table A `data.table` or `data.frame` with at least `cs_name` and `module` columns.
#'        Each row corresponds to a credible set and its module assignment.
#'
#' @return A named list with:
#' \describe{
#'   \item{`master.table`}{The original `module.master.table` with additional columns from
#'   `most.likely.snp.info`, including `most_likely_beta` and `most_likely_se`.}
#'   \item{`joint.lABF.matrix`}{A sparse matrix of joint lABF values, with rows = `cs_name`s and
#'   columns = SNPs. Each `cs_name` receives the lABF profile of its module.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- run_joint_finemapping(
#'   ann.data.obj = ann_data,
#'   module.master.table = master_table
#' )
#'
#' # Updated master table with most likely SNPs
#' head(results$master.table)
#'
#' # Access sparse lABF matrix
#' sparse_mat <- results$joint.lABF.matrix
#' }
#'
#' @seealso \code{\link{joint.finemapping}}
#' 
#' @importFrom data.table as.data.table setnames rbindlist
#' @importFrom Matrix sparseMatrix
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
  
  names(results) <- as.character(names(module.list))
  
  # 1) Extract most likely SNP and joint fine-mapping info, add module reference
  snp_info_list <- list()
  cs_info_list <- list()
  
  for (mod in names(results)) {
    res <- results[[mod]]
    mod_num <- as.numeric(mod)
    
    # Process most.likely.snp.info
    dt_snp <- as.data.table(res$most.likely.snp.info)
    dt_snp[, module := mod_num]
    setnames(dt_snp, "beta", "most_likely_beta")
    setnames(dt_snp, "se", "most_likely_se")
    snp_info_list[[mod]] <- dt_snp
    
    # Process cs.info
    dt_cs <- as.data.table(res$cs.info)
    dt_cs[, module := mod_num]
    cs_info_list[[mod]] <- dt_cs
  }
  
  snp_info_combined <- rbindlist(snp_info_list)
  cs_info_combined <- rbindlist(cs_info_list)

  # 2) Merge most likely SNP info with module.master.table
  module.master.table <- merge(
    module.master.table,
    snp_info_combined,
    by = c("module", "cs_name"),
    all.x = TRUE
  )
  
  # 3) Create sparse matrix: rows = cs_name (from module.master.table), cols = SNPs, values = lABF
  
  # Get mapping of module → cs_names from module.master.table
  cs_name_map <- unique(module.master.table[, .(module, cs_name)])
  
  # Expand cs_info_combined for each cs_name in that module
  cs_info_expanded <- cs_info_combined[cs_name_map, on = "module", allow.cartesian = TRUE] # Now each row of cs_info_combined is duplicated for every cs_name in that module
  
  # Finally build sparse matrix
  row_factors <- factor(cs_info_expanded$cs_name)
  col_factors <- factor(cs_info_expanded$snp)
  
  sparse_mat <- Matrix::sparseMatrix(
    i = as.integer(row_factors),
    j = as.integer(col_factors),
    x = cs_info_expanded$lABF,
    dimnames = list(
      levels(row_factors),  # cs_names
      levels(col_factors)   # SNPs
    )
  )
  
  # Return results
  return(list(
    master.table = module.master.table,
    joint.lABF.matrix = sparse_mat
  ))
}


### WOULD IT BE BETTER TO DO THE EXTRACTION FROM ANNDATA OUTSIDE OF THE LOOP? ALSO FOR RUN_COLOC_PERCOLATION
