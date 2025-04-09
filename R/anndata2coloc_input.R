#' Perform Colocalization Analysis on AnnData Object
#'
#' This function inputs AnnData objec and retrun csv file for coloc nextflow
#' pipeline. It computes the product of the sparse matrix with its transpose to
#' identify shared elements, retrieves trait information,#' and performs
#' colocalization testing on the identified pairs.
#'
#' @param ad An AnnData object containing genetic data.
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
#'
#' @examples
#' \dontrun{
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#'
#' chr22_molQTL_ad <- read_h5ad("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData/HUVEC_chr22_combined_credible_sets.h5ad")
#'
#' chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")
#'
#' chr22_molQTL_ad$obs$phenotype_id <- str_match(chr22_molQTL_ad$obs$cs_name, "chr[0-9]+_([^_]+_[0-9]+|ENSG[0-9]+)")[,2]
#'
#' coloc_guide_table <- anndata2coloc_input(chr22_molQTL_ad)
#'
#' fwrite(coloc_guide_table,file = "/Users/sodbo.sharapov/OneDrive - Htechnopole/00_Sodbo_Projects/anndata/data/coloc_duie_table.csv")
#'}
#'
anndata2coloc_input <- function(ad) {

  matrix_product <- ad$X %*% Matrix::t(ad$X)

  # Set the diagonal to zero (as each vector will have 1s on diagonal)
  diag(matrix_product) <- 0

  shared_elements <- which(matrix_product > 0, arr.ind = TRUE)
  shared_elements_df <- as.data.frame(shared_elements)

  # Map row and column indices to credible set names
  shared_elements_df$row <- rownames(shared_elements)
  shared_elements_df$col <- colnames(matrix_product)[shared_elements_df$col]

  # Keep unique pairs
  shared_elements_unique <- unique(t(apply(shared_elements_df, 1, sort)))
  colnames(shared_elements_unique) <- c("t1", "t2")
  rownames(shared_elements_unique) <- NULL

  # Retrieve trait info
  coloc_combo <- merge(
    shared_elements_unique,
    ad$obs %>% select(-any_of("panel")),
    by.x = "t1",
    by.y = "cs_name"
  )

  coloc_combo <- merge(
    coloc_combo,
    ad$obs %>% select(-any_of("panel")),
    by.x = "t2",
    by.y = "cs_name",
    suffixes = c("_t1", "_t2")
  )

  # Remove pair testing different conditional dataset for the same trait (study_id + phenotype_id)
  coloc_combo <- coloc_combo %>%
    dplyr::filter(study_id_t1 != study_id_t2 | (study_id_t1 == study_id_t2 & phenotype_id_t1 != phenotype_id_t2)) %>%
    dplyr::select(-chr_t2, -contains("start"), -contains("end"))

  coloc_combo <- coloc_combo %>%
    mutate(chr = chr_t1) %>%
    select(
      t2, t1,
      study_id_t1, phenotype_id_t1, any_of("top_pvalue_t1"),
      study_id_t2, phenotype_id_t2, any_of("top_pvalue_t2"),
      chr
    )

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

