#' Create Master Table with Module Assignments
#'
#' @description
#' This function constructs a master table from either a `SingleCellExperiment`
#' or `anndata`-like object, then joins it with a module assignment table.
#' Any rows without an assigned module (i.e., `NA` values in the `module` column)
#' are assigned new unique module IDs starting from the next available integer.
#'
#' @param ad_or_sce Either a `SingleCellExperiment` object or a list mimicking an
#'   `anndata` object (must contain a named element `obs` with cell-level metadata).
#' @param modules A `data.table` containing columns `cs_name` and `module`, representing
#'   precomputed module assignments for a subset of cells or features.
#'
#' @return A `data.table` combining the metadata from `ad_or_sce` with module
#'   labels from `modules`. Any cells/features not assigned to a module will be
#'   given unique module IDs starting from `max(modules$module) + 1`.
#'
#' @details
#' The function determines whether the input object is a `SingleCellExperiment`
#' or an `anndata`-like list, then extracts metadata accordingly. It uses a
#' `data.table` join on the `cs_name` column to merge module labels, and applies
#' the `replace_nas_with_modules()` helper to handle unmatched rows.
#'
#' @examples
#' \dontrun{
#' # Assuming `sce` is a SingleCellExperiment and `modules` a data.table
#' result <- create_master_table(sce, modules)
#'
#' # Or with anndata-like list:
#' result <- create_master_table(list(obs = my_df), modules)
#' }
#'
#' @importFrom SingleCellExperiment colData
#' @import data.table
#' @export
#'
create_master_table <- function(ad_or_sce, modules){
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")
  master.table <- if (is_sce) SingleCellExperiment::colData(ad_or_sce) else ad_or_sce$obs ### anndata obs to create master table
  master.table <- as.data.table(master.table)
  
  master.table_w_modules <- modules[master.table, on = "cs_name"]
  
  ### add module labels single cs modules (NA)
  n_max <- max(modules$module) + 1
  
  ### Full table
  if(any(is.na(master.table_w_modules$module))){
    master.table_w_modules <- replace_nas_with_modules(master.table_w_modules, start_value = n_max)
  }
  return(master.table_w_modules)
}