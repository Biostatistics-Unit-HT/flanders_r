#' Perform Pairwise Colocalization Analysis for a Specific Row
#'
#' This function performs colocalization analysis for a specific pair of traits
#' as defined in a row of the colocalization guide table. It retrieves relevant
#' SNPs within a specified region, imputes missing values, and calculates
#' colocalization probabilities between two datasets.
#'
#' @param coloc_combo_row A row from the colocalization guide table containing
#' information on trait pairs to be tested.
#'
#' @param ad_by_chr A list of AnnData object, split by chromosome. Each element
#' of the list should have name like chr1, chr 21 and so on
#'
#' @return A list with colocalization results, including a summary with the top
#' SNPs and trait information.
#' @export
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' coloc_input <- anndata2coloc_input(ad)
#' coloc_combo_row <- coloc_input[1,]
#' # Example usage of anndata2coloc_row
#' result <- anndata2coloc_row(coloc_combo_row)
#' }

anndata2coloc_row <- function(coloc_combo_row, ad_by_chr) {

  chr_name <- as.character(ad$obs[coloc_combo_row[1,]$t1,"chr"])

  ad_sm <- ad_by_chr[[chr_name]]

  loc_start <- min(ad_sm$obs[coloc_combo_row$t1,"start"],ad_sm$obs[coloc_combo_row$t2,"start"])
  loc_start <- loc_start - 1e5

  loc_end <- max(ad_sm$obs[coloc_combo_row$t1,"end"],ad_sm$obs[coloc_combo_row$t2,"end"])
  loc_end <- loc_end + 1e5

  # Load-in precomputed lABF
  snps_for_coloc <- ad_sm$var %>%
    filter(
      pos >= loc_start,
      pos <= loc_end
    ) %>% rownames

  dataset <- ad_sm[c(coloc_combo_row$t1,coloc_combo_row$t2), snps_for_coloc]$X

  dataset1 <- dataset[1,]
  names(dataset1) <- snps_for_coloc
  #dataset1 <- dataset1[which(dataset1 != 0)]

  dataset2 <- dataset[2,]
  names(dataset2) <- snps_for_coloc
  #dataset2 <- dataset2[which(dataset2 != 0)]

  #dataset1 <- dataset1[which(names(dataset1) %in% snps_for_coloc)]
  #dataset2 <- dataset2[which(names(dataset2) %in% snps_for_coloc)]

  dataset1 <- impute.labf(
    snp = snps_for_coloc,
    cred.set = names(dataset1)[dataset1 != 0],
    cred.set.labf = dataset1[dataset1 != 0],
    imputed.labf = ad[coloc_combo_row$t1,]$obs[1, "min_res_labf"]
  )

  dataset2 <- impute.labf(
    snp = snps_for_coloc,
    cred.set = names(dataset2)[dataset2 != 0],
    cred.set.labf = dataset2[dataset2 != 0],
    imputed.labf = ad[coloc_combo_row$t2,]$obs[1,"min_res_labf"]
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
      t1_study_id = coloc_combo_row$t1_study_id,
      t1 = t1,
      t2_study_id = coloc_combo_row$t2_study_id,
      t2 = t2,
      hit1 = cojo_snp1,
      hit2 = cojo_snp2
    )

  icoloc.res

}
