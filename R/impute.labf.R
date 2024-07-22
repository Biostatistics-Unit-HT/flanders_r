#' Impute lABF Values for SNPs
#'
#' This function imputes lABF values for SNPs. For SNPs in the credible set, it assigns their corresponding lABF values.
#' For SNPs not in the credible set, it assigns an imputed lABF value.
#'
#' @param snp A character vector with names of SNPs.
#' @param cred.set A character vector of SNPs which are in the credible set.
#' @param cred.set.labf A numeric vector of lABF values for SNPs in the credible set.
#' @param imputed.labf A numeric value which will be used for imputation of lABF values for SNPs not in the credible set.
#'
#' @return A data frame with SNP names and their corresponding lABF values, either observed or imputed.
#' @export
#'
#' @examples
#' \dontrun{
#' snps <- c("rs1", "rs2", "rs3", "rs4")
#' credible_set <- c("rs1", "rs3")
#' credible_set_labf <- c(10.5, 12.3)
#' imputed_labf_value <- 5.0
#' impute.labf(snps, credible_set, credible_set_labf, imputed_labf_value)
#' }

impute.labf <- function(snp, cred.set, cred.set.labf, imputed.labf) {
  
  # Create a data frame to store the SNPs and their lABF values
  df_imp <- data.frame(
    snp = snp,
    lABF = NA,
    row.names = snp
  )
  
  # Assign lABF values for SNPs in the credible set
  df_imp[match(cred.set, snp), "lABF"] <- cred.set.labf
  
  # Identify SNPs with NA lABF values
  markers_to_impute <- is.na(df_imp$lABF)
  
  # Impute lABF values for SNPs not in the credible set
  df_imp$lABF[markers_to_impute] <- imputed.labf
  
  return(df_imp)
}
