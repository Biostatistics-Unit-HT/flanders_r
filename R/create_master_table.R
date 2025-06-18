#' Create Master Table with Module Assignments
#'
#' @description
#' This function constructs a master table from either a `SingleCellExperiment`
#' or `anndata`-like object, then joins it with module assignment tables.
#' The function supports either a single set of simple modules or a combination
#' of post-coloc percolation and fully connected modules. Any rows without an
#' assigned module (i.e., `NA` values in the `module` column) are assigned new
#' unique module IDs starting from the next available integer.
#'
#' @param ad_or_sce Either a `SingleCellExperiment`  or an `anndata`-like object
#'(must contain a named element `obs` with credible set-level metadata).
#' @param simple_modules A data.table or data.frame containing simple module assignments. 
#'        Must include a `cs_name` and `module` column. Optional.
#' @param partially_connected_modules A data.table or data.frame with post-coloc percolation
#'        module assignments. Used when `simple_modules` is NULL. Must include `cs_name` and 
#'        string-based `module` values. Optional.
#' @param fully_connected_modules A data.table or data.frame with numeric module assignments 
#'        for fully connected modules. Used when `simple_modules` is NULL. Optional.
#'
#' @return A `data.table` combining credible set metadata with numeric module assignments. 
#'         NA module values are replaced with unique numeric IDs.
#'
#' @details If `simple_modules` is provided, it is used directly. Otherwise, both 
#' `partially_connected_modules` and `fully_connected_modules` must be provided. The function 
#' will remap string module identifiers to numeric values to avoid overlap.
#'
#' @examples
#' \dontrun{
#' # Using SingleCellExperiment
#' master_table <- create_master_table(
#'   ad_or_sce = sce,
#'   simple_modules = module_df
#' )
#'
#' # Using anndata (converted to R list via reticulate)
#' master_table <- create_master_table(
#'   ad_or_sce = ad,
#'   partially_connected_modules = pc_modules,
#'   fully_connected_modules = fc_modules
#' )
#' }
#'
#' @importFrom SingleCellExperiment colData
#' @import data.table
#' @export
#'
create_master_table <- function(
    ad_or_sce,
    simple_modules = NULL,
    partially_connected_modules = NULL,
    fully_connected_modules = NULL
) {
  
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")
  master.table <- if (is_sce) SingleCellExperiment::colData(ad_or_sce) else ad_or_sce$obs ### anndata obs to create master table
  master.table <- as.data.table(master.table)
  
  # Handle module inputs
  if (!is.null(simple_modules)) {
    combined_modules <- as.data.table(simple_modules)
  } else if (!is.null(post_cp_modules) && !is.null(fully_connected_modules)) {
    post_cp_modules <- as.data.table(post_cp_modules)
    fully_connected_modules <- as.data.table(fully_connected_modules)
    
    # Ensure fully_connected_modules$module is numeric
    fully_connected_modules[, module := as.integer(module)]
    
    # Get max value from fully_connected_modules to avoid collisions
    offset <- max(fully_connected_modules$module, na.rm = TRUE) + 1
    
    # Map string modules in partially_connected_modules to numeric values starting from offset
    unique_modules <- unique(post_cp_modules$module)
    module_mapping <- setNames(seq(offset, offset + length(unique_modules) - 1), unique_modules)
    post_cp_modules[, module := module_mapping[module]]
    
    # Combine both
    combined_modules <- rbindlist(list(fully_connected_modules, post_cp_modules), use.names = TRUE, fill = TRUE)
  } else {
    stop("Provide either `modules`, or both `modules1` and `modules2`.")
  }
  
  # Step 2: Merge with metadata
  master.table_w_modules <- combined_modules[master.table, on = "cs_name"]
  
  # Step 3: Handle NA modules — assign unique numeric ID per NA row
  next_module <- max(combined_modules$module, na.rm = TRUE) + 1
  
  if (any(is.na(master.table_w_modules$module))) {
    master.table_w_modules <- replace_nas_with_modules(
      master.table_w_modules,
      start_value = next_module
    )
  }
  
  return(master.table_w_modules)
}
