#' Perform Colocalization Analysis on AnnData Object
#'
#' This function inputs AnnData objec and retrun csv file for coloc nextflow
#' pipeline. It computes the product of the sparse matrix with its transpose to
#' identify shared elements, retrieves trait information,#' and performs
#' colocalization testing on the identified pairs.
#'
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing genetic data.
#'
#' @return A data.frame with colnames of t2, t1, t1_study_id,t1_phenotype_id,
#' t1_top_pvalue, t2_study_id, t2_phenotype_id, t2_top_pvalue. Each row
#' represents a pariwise colocalization test which is needed to be completed.
#' This aata.frame, written on hard drive as CSV file is an input for
#' anndata2coloc function.
#'
#' @export
#' @import data.table
#' @import dplyr
#' @import anndata
#' @import coloc
#' @import Matrix
#' @import zellkonverter
#' @import SingleCellExperiment
#' @import scRNAseq
#'
#' @title Convert AnnData or SingleCellExperiment to Coloc Input
#' 
#' @description This script demonstrates how to process credible sets data 
#' stored in an AnnData or SingleCellExperiment object to generate input for 
#' coloc analysis. It extracts relevant metadata, such as study IDs and 
#' phenotype IDs, and formats the data into a table suitable for coloc.
#' 
#' @details 
#' The script assumes that the input data is stored in an AnnData object 
#' (e.g., `.h5ad` file) or a SingleCellExperiment object. It extracts the 
#' `study_id` and `phenotype_id` from the `cs_name` field in the `obs` 
#' (observations) metadata. The resulting table is saved as a CSV file for 
#' downstream coloc analysis.
#' 
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing 
#' credible sets data. The `obs` metadata should include a `cs_name` field.
#' 
#' @return A data frame formatted for coloc analysis, saved as a CSV file.
#' 
#' @examples
#' \dontrun{
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#'
#' # Load AnnData object
#' chr22_molQTL_ad <- read_h5ad("/path/to/HUVEC_chr22_combined_credible_sets.h5ad")
#'
#' # Extract study and phenotype IDs
#' chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")
#' chr22_molQTL_ad$obs$phenotype_id <- str_match(
#'   chr22_molQTL_ad$obs$cs_name,
#'   "chr[0-9]+_([^_]+_[0-9]+|ENSG[0-9]+)"
#' )[,2]
#'
#' # Generate coloc input table
#' coloc_guide_table <- anndata2coloc_input(chr22_molQTL_ad)
#'
#' # Save the table to a CSV file
#' fwrite(coloc_guide_table, file = "/path/to/coloc_guide_table.csv")
#' }
#' 
#' @seealso 
#' \code{\link[anndata]{read_h5ad}} for reading AnnData files.
#' \code{\link[data.table]{fwrite}} for writing data to CSV files.

anndata2coloc_input <- function(ad_or_sce) {
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")

  if (is_sce) {
    # ---- Case: SCE object from zellkonverter ----
    message("Processing SingleCellExperiment object...")

    X <- as(Matrix::t(assay(ad_or_sce, "X")), "dgCMatrix")
    cs_names <- rownames(colData(ad_or_sce))
    obs_df <- as.data.frame(colData(ad_or_sce))
  } else {
    # ---- Case: AnnData object from reticulate ----
    message("Processing AnnData object from Python...")

    X <- ad_or_sce$X
    cs_names <- rownames(ad_or_sce$obs)
    obs_df <- as.data.frame(ad_or_sce$obs)
  }

  # Matrix multiplication to detect shared SNPs
  matrix_product <- X %*% Matrix::t(X)
  diag(matrix_product) <- 0

  # Find non-zero overlaps (shared SNPs)
  shared_elements <- which(matrix_product > 0, arr.ind = TRUE)
  shared_elements_df <- as.data.frame(shared_elements)
  shared_elements_df$row <- cs_names[shared_elements_df$row]
  shared_elements_df$col <- cs_names[shared_elements_df$col]

  # Get unique and sorted credible set pairs
  shared_elements_unique <- unique(t(apply(shared_elements_df[, c("row", "col")], 1, sort)))
  colnames(shared_elements_unique) <- c("t1", "t2")
  rownames(shared_elements_unique) <- NULL

  # Merge trait info
  coloc_combo <- merge(
    shared_elements_unique,
    obs_df %>% dplyr::select(-any_of("panel")),
    by.x = "t1",
    by.y = "cs_name"
  )

  coloc_combo <- merge(
    coloc_combo,
    obs_df %>% dplyr::select(-any_of("panel")),
    by.x = "t2",
    by.y = "cs_name",
    suffixes = c("_t1", "_t2")
  )

  # Remove same study_id and phenotype_id duplicates
  coloc_combo <- coloc_combo %>%
    dplyr::filter(
      study_id_t1 != study_id_t2 |
        (study_id_t1 == study_id_t2 & phenotype_id_t1 != phenotype_id_t2)
    ) %>%
    dplyr::select(-chr_t2, -contains("start"), -contains("end"))

  coloc_combo <- coloc_combo %>%
    mutate(chr = chr_t1) %>%
    select(
      t2, t1,
      study_id_t1, phenotype_id_t1, any_of("top_pvalue_t1"),
      study_id_t2, phenotype_id_t2, any_of("top_pvalue_t2"),
      chr
    )

  # Rename columns based on available info
  if (ncol(coloc_combo) == 9) {
    colnames(coloc_combo) <- c(
      "t2", "t1",
      "t1_study_id", "t1_phenotype_id", "t1_top_pvalue",
      "t2_study_id", "t2_phenotype_id", "t2_top_pvalue",
      "chr"
    )
  } else if (ncol(coloc_combo) == 7) {
    colnames(coloc_combo) <- c(
      "t2", "t1",
      "t1_study_id", "t1_phenotype_id",
      "t2_study_id", "t2_phenotype_id",
      "chr"
    )
  } else {
    stop("Unexpected number of columns in coloc_combo")
  }

  return(coloc_combo)
}