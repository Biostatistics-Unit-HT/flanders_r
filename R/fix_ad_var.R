#' Fix Variant Information in AnnData Object
#'
#' This function modifies the `var` component of an AnnData object by adding a
#' `snp` column with the row names of the `var` dataframe. It also extracts
#' chromosome (`chr`) and position (`pos`) information from the `snp`
#' identifiers and assigns them to new columns in `var`.
#'
#' @param ad An AnnData object containing genetic variant data in its `var`
#' component.
#'
#' @return The modified AnnData object with updated `snp`, `chr`, and `pos`
#' columns in the `var` component.
#' @export
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(anndata)
#' library(dplyr)
#'
#' # Load AnnData object
#' ad <- read_h5ad("path/to/anndata_file.h5ad")
#'
#' # Fix variant information in the AnnData object
#' ad <- fix_ad_var(ad)
#'
#' # Check the updated var component
#' head(ad$var)
#' }
fix_ad_var <- function(ad) {

  ad$var$snp <- rownames(ad$var)

  ad$var <- ad$var %>%
    mutate(
      chr = as.character(sub("(:[0-9]+:.*)", "", snp)),
      pos = as.numeric(sub("chr[0-9]+:([0-9]+):.*", "\\1", snp))
    )

  ad

}
