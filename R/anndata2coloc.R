#' Perform Colocalization Analysis on AnnData Object
#'
#' This function performs colocalization analysis on an AnnData object. It computes the product of
#' the sparse matrix with its transpose to identify shared elements, retrieves trait information,
#' and performs colocalization testing on the identified pairs.
#'
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing genetic data.
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
anndata2coloc <- function(ad_or_sce, coloc_input) {
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")

  # Ensure chromosome columns are character
  if (is_sce) {
    rowData(ad_or_sce)$chr <- as.character(rowData(ad_or_sce)$chr)
    colData(ad_or_sce)$chr <- as.character(colData(ad_or_sce)$chr)
  } else {
    ad_or_sce$obs$chr <- as.character(ad_or_sce$obs$chr)
    ad_or_sce$var$chr <- as.character(ad_or_sce$var$chr)
  }
    
  # Extract chromosome information:
  # For SCE: credible set info is in colData,
  # and the credible set names are in the rownames of colData (which should equal colnames(ad_or_sce)).
  # For AnnData (reticulate): assume the same structure in ad_or_sce$obs.
  obs_chr    <- if (is_sce) colData(ad_or_sce)$chr else ad_or_sce$obs$chr
  obs_csname <- if (is_sce) rownames(colData(ad_or_sce)) else rownames(ad_or_sce$obs)
  
  # Determine which chromosomes the input credible sets (coloc_input$t1) belong to
  coloc_chr  <- obs_chr[match(coloc_input$t1, obs_csname)]
  list_of_chr <- unique(coloc_chr)
  
  message("There are ", length(list_of_chr), " chromosomes: ", paste(list_of_chr, collapse = ", "))
  
  # Subset the object by chromosome.
  # For SCE: rows = SNPs (from rowData, also rownames(ad_or_sce)) and columns = credible sets (from colData, also colnames(ad_or_sce))
  ad_by_chr <- lapply(list_of_chr, function(chr) {
    if (is_sce) {
      snp_rows  <- which(rowData(ad_or_sce)$chr == chr)   # SNPs (features)
      cred_cols <- which(colData(ad_or_sce)$chr == chr)     # Credible sets (cells)
      ad_or_sce[snp_rows, cred_cols]
    } else {
      ad_or_sce[ad_or_sce$obs$chr == chr, ad_or_sce$var$chr == chr]
    }
  })
  names(ad_by_chr) <- as.character(list_of_chr)
  
  # Split coloc_input into one-row pieces
  coloc_combo_ls <- split(coloc_input, seq_len(nrow(coloc_input)))
  
  # For each input row, perform colocalization analysis
  coloc.full <- lapply(coloc_combo_ls, anndata2coloc_row, ad_by_chr = ad_by_chr)
  
  # Combine the individual summaries into one data.frame
  only_summary_df <- data.table::rbindlist(
    lapply(coloc.full, function(x) x$summary),
    fill = TRUE
  )
  
  return(as.data.frame(only_summary_df))
}
