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
#' chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")
#'
#' chr22_molQTL_ad$obs$phenotype_id <- str_match(chr22_molQTL_ad$obs$cs_name, "chr[0-9]+_([^_]+_[0-9]+|ENSG[0-9]+)")[,2]
#'
#' coloc_guide_table <_ anndata2coloc_input(chr22_molQTL_ad)
#'
#' coloc_res <- anndata2coloc(chr22_molQTL_ad,coloc_guide_table)
#'
#'
anndata2coloc <- function(ad, coloc_input) {

  ###### COLOCALIZATION ######

  # Split the dataframe into a list of rows
  coloc_combo_ls <- split(coloc_input, seq(nrow(coloc_input)))

  # Perform coloc!
  coloc.full <- lapply(coloc_combo_ls, function(coloc_combo_row) {

    # Load-in precomputed lABF
    dataset1 <- ad$X[coloc_combo_row$t1,]
    dataset1 <- dataset1[which(dataset1 != 0)]

    dataset2 <- ad$X[coloc_combo_row$t2,]
    dataset2 <- dataset2[which(dataset2 != 0)]

    dataset1 <- impute.labf(
      snp = union(names(dataset1), names(dataset2)),
      cred.set = names(dataset1),
      cred.set.labf = dataset1,
      imputed.labf = ad$obs[coloc_combo_row$t1, "min_res_labf"]
    )

    dataset2 <- impute.labf(
      snp = union(names(dataset2), names(dataset2)),
      cred.set = names(dataset2),
      cred.set.labf = dataset2,
      imputed.labf = ad$obs[coloc_combo_row$t2, "min_res_labf"]
    )

    icoloc.res <- icoloc.ht(
      dataset1,
      dataset2
    )

    # Retrieve important info from file name
    t1 <- ifelse(is.na(coloc_combo_row$t1_phenotype_id), coloc_combo_row$t1_study_id, paste0(coloc_combo_row$t1_study_id, "_", coloc_combo_row$t1_phenotype_id))
    cojo_snp1 <- gsub(paste0(t1, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t1)

    t2 <- ifelse(is.na(coloc_combo_row$t2_phenotype_id), coloc_combo_row$t2_study_id, paste0(coloc_combo_row$t2_study_id, "_", coloc_combo_row$t2_phenotype_id))
    cojo_snp2 <- gsub(paste0(t2, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t2)

    # Add top cojo SNPs and traits
    icoloc.res$summary <- icoloc.res$summary %>%
      mutate(
        t1_study_id = coloc_combo_row$t1_sty_id,
        t1 = t1,
        t2_study_id = coloc_combo_row$t2_study_id,
        t2 = t2,
        hit1 = cojo_snp1,
        hit2 = cojo_snp2
      )

    icoloc.res

  })

# Store ALL the summary output in a data frame, adding tested traits column and SAVE
  only_summary_df <- as.data.frame(rbindlist(lapply(coloc.full, function(x) { x$summary })))

  return(only_summary_df)

}
