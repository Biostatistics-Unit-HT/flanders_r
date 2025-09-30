#' Fix Variant Information in AnnData or SingleCellExperiment Object
#'
#' This function standardizes variant metadata by ensuring that a `snp` column 
#' exists and extracting chromosome (`chr`) and position (`pos`) information 
#' from SNP identifiers. Works with both Python AnnData objects (via \pkg{anndata}) 
#' and R \code{SingleCellExperiment} objects (via \pkg{zellkonverter}).
#'
#' SNP identifiers are expected in the format: \code{"chr{num}:{pos}:EA:RA"}.
#'
#' @param ad An object of class \code{AnnData} (Python, via reticulate) 
#'   or \code{SingleCellExperiment} (R, via Bioconductor).
#'
#' @return The modified object, with updated `snp`, `chr`, and `pos` columns 
#'   in the variant metadata:
#'   \itemize{
#'     \item For \code{AnnData}: in \code{ad$var}.
#'     \item For \code{SingleCellExperiment}: in \code{rowData(ad)}.
#'   }
#'
#' @export
#' @import dplyr
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#' \dontrun{
#' library(anndata)
#' library(zellkonverter)
#' library(SingleCellExperiment)
#'
#' # Example 1: AnnData object (Python via reticulate)
#' ad <- read_h5ad("path/to/anndata_file.h5ad")
#' ad <- fix_ad_var(ad)
#' head(ad$var)
#'
#' # Example 2: SingleCellExperiment object
#' sce <- readH5AD("path/to/anndata_file.h5ad", reader = "R")
#' sce <- fix_ad_var(sce)
#' head(rowData(sce))
#' }
fix_ad_var <- function(ad) {
  if (inherits(ad, "SingleCellExperiment")) {
    # Extract rowData and ensure it's a DataFrame
    var <- as.data.frame(rowData(ad))
    
    # Create snp column from rownames if missing
    if (!"snp" %in% colnames(var)) {
      var$snp <- rownames(ad)
    }
    
    # Add chr and pos columns
    var <- var %>%
      mutate(
        chr = sub(":.*", "", snp),
        pos = as.numeric(sub(".*:(\\d+):.*", "\\1", snp))
      )
    
    # Write back into rowData
    rowData(ad) <- S4Vectors::DataFrame(var)
    
  } else if (inherits(ad, "AnnData")) {
    # Add snp column if missing
    if (!"snp" %in% colnames(ad$var)) {
      ad$var$snp <- rownames(ad$var)
    }
    
    ad$var <- ad$var %>%
      mutate(
        chr = sub(":.*", "", snp),
        pos = as.numeric(sub(".*:(\\d+):.*", "\\1", snp))
      )
  } else {
    stop("Unsupported object type: must be AnnData or SingleCellExperiment")
  }
  
  ad
}
