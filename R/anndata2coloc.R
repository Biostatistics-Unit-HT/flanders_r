#' Perform Colocalization Analysis on AnnData Object
#'
#' This function performs colocalization analysis on an AnnData object. It computes the product of
#' the sparse matrix with its transpose to identify shared elements, retrieves trait information,
#' and performs colocalization testing on the identified pairs.
#'
#' @param ad An AnnData object containing genetic data.
#' @param coloc_input a path to a file with schema of pairwises test to run
#'
#' @return A list of colocalization results for each pair of traits.
#' @export
#' @import data.table
#' @import dplyr
#' @import anndata
#' @import coloc
#' @import Matrix
#'
#' @examples
#' \dontrun{
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#'
#' chr22_molQTL_ad <- read_h5ad("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData/HUVEC_chr22_combined_credible_sets.h5ad")
#'
#' coloc_guide_table <- anndata2coloc_input(chr22_molQTL_ad)
#'
#' coloc_res <- anndata2coloc(chr22_molQTL_ad,coloc_guide_table)
#'
#'}
#'
anndata2coloc <- function(ad, coloc_input) {

  list_of_chr <- unique(ad$obs[coloc_input$t1,"chr"])

  print(paste("There are",list_of_chr))

  ad_by_chr <- lapply(list_of_chr,function(chr) {ad[ad$obs$chr == chr, ad$var$chr == chr]})

  names(ad_by_chr) <- as.character(list_of_chr)

  ###### COLOCALIZATION ######

  # Split the dataframe into a list of rows
  coloc_combo_ls <- split(coloc_input, seq(nrow(coloc_input)))

  # Perform coloc!
  coloc.full <- lapply(coloc_combo_ls, anndata2coloc_row,ad_by_chr = ad_by_chr)

# Store ALL the summary output in a data frame, adding tested traits column and SAVE
  only_summary_df <- as.data.frame(rbindlist(lapply(coloc.full, function(x) { x$summary })))

  return(only_summary_df)

}
