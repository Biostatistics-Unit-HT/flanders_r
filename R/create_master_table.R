#' Create Master Table with Module Assignments
#'
#' @description
#' This function constructs a master table from either a `SingleCellExperiment`
#' or `anndata`-like object, then joins it with module assignment tables.
#' It supports either a single set of simple modules or a combination
#' of refined (post-coloc percolation or hypercoloc) and fully connected modules.
#' When the output list of `run_coloc_percolation()` or `community2hypercoloc()`
#' is provided, the `post_cp_tmp_module` table is extracted from each element.
#' In addition, any rows without an assigned module 
#' (i.e., `NA` values in the `module` column) are assigned new
#' unique module IDs starting from the next available integer.
#'
#' @param ad_or_sce Either a `SingleCellExperiment`  or an `anndata`-like object
#' (must contain a named element `obs` with credible set-level metadata).
#' @param simple_modules A data.table or data.frame containing simple module assignments. 
#'        Must include a `cs_name` and `module` column. Optional.
#' @param post_cp_list A **list** object where each element contains `post_cp_tmp_module`
#'        (e.g., output from `run_coloc_percolation()`). Optional.
#' @param fully_connected_modules A data.table or data.frame with numeric module assignments 
#'        for fully connected modules. Optional.
#'
#' @return A `data.table` combining credible set metadata with numeric module assignments. 
#'         NA module values are replaced with unique numeric IDs.
#'
#' @details If `simple_modules` is provided, it is used directly. Otherwise, both 
#' `post_cp_list` and `fully_connected_modules` must be provided. The function 
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
#'   post_cp_list = coloc_percolation_groups,
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
    post_cp_list = NULL,
    fully_connected_modules = NULL
) {
  
# Set up master table  
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")
  master.table <- if (is_sce) SingleCellExperiment::colData(ad_or_sce) else ad_or_sce$obs
  master.table <- as.data.table(master.table)
  
  # Handle module inputs
  if (!is.null(simple_modules)) {
    
    combined_modules <- as.data.table(simple_modules)
    combined_modules[, module_type := "simple"]
    
  } else if (!is.null(post_cp_list) && !is.null(fully_connected_modules)) {
    
    # Extract and combine all post-coloc modules
    post_cp_modules <- rbindlist(
      lapply(post_cp_list, function(x) x$post_cp_tmp_module),
      use.names = TRUE)
    post_cp_modules <- as.data.table(post_cp_modules)
    
    fully_connected_modules <- as.data.table(fully_connected_modules)
    
    # Ensure fully_connected_modules$module is numeric
    fully_connected_modules[, module := as.integer(module)]
    
    # Add type of module
    fully_connected_modules[, module_type := "fully_connected"]
    
    # Get max value from fully_connected_modules to avoid collisions
    offset <- max(fully_connected_modules$module, na.rm = TRUE) + 1
    
    # Map string modules in partially_connected_modules to numeric values starting from offset
    unique_modules <- unique(post_cp_modules$module)
    module_mapping <- setNames(seq(offset, offset + length(unique_modules) - 1), unique_modules)
    post_cp_modules[, module := module_mapping[module]]
    
    # Add type of module
    post_cp_modules[, module_type := "refined"]
    
    # Combine both
    combined_modules <- rbindlist(list(fully_connected_modules, post_cp_modules), use.names = TRUE, fill = TRUE)
  } else {
    stop("Provide either `simple_modules`, or both `post_cp_list` and `fully_connected_modules`.")
  }
  
  # Step 2: Merge with metadata
  master.table_w_modules <- combined_modules[master.table, on = "cs_name"]
  master.table_w_modules$module_type <- ifelse(is.na(master.table_w_modules$module_type), "fully_connected", master.table_w_modules$module_type)
  
  # Step 3: Handle NA modules — assign unique numeric ID per NA row
  next_module <- max(combined_modules$module, na.rm = TRUE) + 1
  
  if (any(is.na(master.table_w_modules$module))) {
    master.table_w_modules <- flanders:::replace_nas_with_modules(
      master.table_w_modules,
      start_value = next_module
    )
  }
  
  return(master.table_w_modules)
}
