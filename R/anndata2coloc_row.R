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
#' result <- anndata2coloc_row(coloc_combo_row, ad_by_chr)
#' }

anndata2coloc_row <- function(coloc_combo_row, ad_by_chr) {
  # Determine if we are working with a SingleCellExperiment (SCE)
  is_sce <- inherits(ad_by_chr[[1]], "SingleCellExperiment")
  chr_name <- as.character(coloc_combo_row$chr[1])
  
  # Obtain the chromosome-specific subset
  ad_sm <- ad_by_chr[[chr_name]]
  
  # For AnnData: use ad_sm$obs and ad_sm$var; for SCE: use colData and rowData.
  ad_obs <- if (is_sce) colData(ad_sm) else ad_sm$obs
  ad_var <- if (is_sce) rowData(ad_sm) else ad_sm$var
  
  t1 <- coloc_combo_row$t1
  t2 <- coloc_combo_row$t2
  
  # Extract start and end positions for each credible set
  start1 <- ad_obs[t1, "start"]
  start2 <- ad_obs[t2, "start"]
  end1   <- ad_obs[t1, "end"]
  end2   <- ad_obs[t2, "end"]
  
  loc_start <- min(start1, start2) - 1e5
  loc_end   <- max(end1, end2) + 1e5
  
  # Identify SNPs in the interval using SNP metadata (ad_var)
  snp_pos   <- ad_var$pos
  snp_names <- rownames(ad_var)
  snps_for_coloc <- snp_names[snp_pos >= loc_start & snp_pos <= loc_end]
  
  # Subset the matrix containing precomputed lABF values
  if (is_sce) {
    # For SCE: 
    # - Features (SNPs) are stored in rowData (and rownames(ad_sm) are SNP names),
    # - Credible sets are stored in colData (and rownames(colData(ad_sm)) are credible set names).
    rows <- match(snps_for_coloc, rownames(rowData(ad_sm)))
    cols <- match(c(t1, t2), rownames(colData(ad_sm)))
    
    if (any(is.na(rows))) {
      warning("Some SNPs in the interval are missing in rowData; these will be omitted.")
      valid <- !is.na(rows)
      snps_for_coloc <- snps_for_coloc[valid]
      rows <- rows[valid]
    }
    if (any(is.na(cols))) {
      stop("Credible set names ", paste(c(t1, t2)[is.na(cols)], collapse = ", "),
           " not found in colData.")
    }
    
    mat <- assay(ad_sm, "X")[rows, cols, drop = FALSE]
    
  } else {
    # For AnnData:
    # The subset operation returns a matrix with dimensions transposed relative to SCE.
    # So we transpose it here to get rows = credible sets and columns = SNPs.
    mat <- t(ad_sm[c(t1, t2), snps_for_coloc]$X)
  }
  
  # Extract lABF vectors for the two credible sets
  dataset1 <- mat[, 1]
  dataset2 <- mat[, 2]
  
  # Ensure that each dataset gets SNP names as names
  names(dataset1) <- snps_for_coloc
  names(dataset2) <- snps_for_coloc
  
  # Retrieve imputation reference (minimum residual lABF) for each credible set
  min_res_labf_1 <- ad_obs[t1, "min_res_labf"]
  min_res_labf_2 <- ad_obs[t2, "min_res_labf"]
  
  # Impute lABF values for SNPs not present in the credible set
  dataset1 <- impute.labf(
    snp = snps_for_coloc,
    cred.set = names(dataset1)[dataset1 != 0],
    cred.set.labf = dataset1[dataset1 != 0],
    imputed.labf = min_res_labf_1
  )
  
  dataset2 <- impute.labf(
    snp = snps_for_coloc,
    cred.set = names(dataset2)[dataset2 != 0],
    cred.set.labf = dataset2[dataset2 != 0],
    imputed.labf = min_res_labf_2
  )
  
  # Run colocalization analysis (using your predefined icoloc.ht function)
  icoloc.res <- icoloc.ht(dataset1, dataset2)
  
  # Format the trait labels and extract a key SNP name from the credible set filename.
  t1_label <- ifelse(
    is.na(coloc_combo_row$t1_phenotype_id),
    coloc_combo_row$t1_study_id,
    coloc_combo_row$t1_phenotype_id
  )
  
  t2_label <- ifelse(
    is.na(coloc_combo_row$t2_phenotype_id),
    coloc_combo_row$t2_study_id,
    coloc_combo_row$t2_phenotype_id
  )
  
  cojo_snp1 <- gsub(paste0(t1_label, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t1)
  cojo_snp2 <- gsub(paste0(t2_label, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t2)
  
  icoloc.res$summary <- icoloc.res$summary %>%
    dplyr::mutate(
      t1_study_id = coloc_combo_row$t1_study_id,
      t1 = t1_label,
      t2_study_id = coloc_combo_row$t2_study_id,
      t2 = t2_label,
      hit1 = cojo_snp1,
      hit2 = cojo_snp2
    )
  
  return(icoloc.res)
}
