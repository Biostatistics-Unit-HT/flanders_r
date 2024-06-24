# Sodbo Sharapov, Human Technopole (c) 2024

#' Filters SuSiE output
#'
#' @description
#' This function filters the output from SuSiE and outputs Susie class object.
#' It return the object of class "susie" with filtered:
#' susie_output$sets$cs
#' susie_output$sets$purity
#' susie_output$sets$cs_index
#' susie_output$sets$coverage
#'
#' @param susie_output A SuSiE_R output (e.g. susieR::fitted_rss()).
#' @param cs_lbf_thr A numeric value representing the credible set log_e_BF threshold for filtering credible sets (default is 2).
#' @param pval A numeric vector representing P-values of the SNPs in the locus, tested in the susie_output.
#' @param cs_lbf_thr A numeric value representing the credible set logBF threshold for filtering credible sets (default is 2).
#' @param signal_pval_threshold A numeric value representing top SNP P-value threshold for filtering credible sets (default is 1e-6).
#' @param purity_mean_r2_threshold A numeric value representing the threshold for purity mean r2 QC metrics for filtering credible sets (default is 0.5).
#' @param purity_min_r2_threshold A numeric value representing the threshold for purity min r2 QC metrics for filtering credible sets (default is 0.5).
#' @param verbose A boolean value. If TRUE, the log information will be printed.
#'
#' @return A Susie object containing filtered fine-mapped credible sets.
#'
#' @details
#' This function takes the output from SuSiE, filters credible sets by purity metrics and CS lbf,
#' and returns a Susie object.
#'
#' @examples
#' \dontrun{
#' susie_filtered <- susie.cs.ht(
#'     susie_output = susie_output,
#'     pval = pval,
#'     cs_lbf_thr = 2,
#'     signal_pval_threshold = 1,
#'     purity_mean_r2_threshold = 0,
#'     purity_min_r2_threshold = 0,
#'     verbose = TRUE
#' )
#' }
#' @export
susie.cs.ht <- function(
    susie_output = NULL,
    pval = NULL,
    cs_lbf_thr = 2,
    signal_pval_threshold = 1,
    purity_mean_r2_threshold = 0.5,
    purity_min_r2_threshold = 0.5,
    verbose = FALSE
) {

  # Head of summary.susie()$cs:
  # cs cs_log10bf cs_avg_r2 cs_min_r2 variable
  # 2  22.803334 1.0000000 1.0000000 371
  # 3  13.595152 1.0000000 1.0000000 255
  # 8   4.890725 1.0000000 1.0000000 492
  # 2.738234 0.9994622 0.9994622 200,202

  cs_info <- susieR::summary.susie(susie_output)$cs

  # Get names of the CS that have lbf above threshold
  good_lbf_cs <- cs_info$cs[cs_info$cs_log10bf >= log10(exp(cs_lbf_thr))]

  # Get names of the CS that have mean purity above threshold
  good_purity_mean_cs <-  cs_info$cs[cs_info$cs_avg_r2 >= purity_mean_r2_threshold]

  # Get names of the CS that have min purity above threshold
  good_purity_min_cs <-  cs_info$cs[cs_info$cs_min_r2 >= purity_min_r2_threshold]

  if(!is.null(pval)) {
    # Get vector of indexes of SNPs per each CS
    min_pval_per_cs <- lapply(
      strsplit(cs_info$variable,split = ","),
      as.numeric)

    # Get vector of P-values of top SNP per each CS
    min_pval_per_cs <- sapply(min_pval_per_cs,function(x) min(pval[x]))

    # Get names of the CS that have index SNP with P-value below threshold
    good_pval_cs <- cs_info$cs[min_pval_per_cs <= signal_pval_threshold]

  } else {

    good_pval_cs <- cs_info$cs

  }

  # Get intersect of CS that satisfy all the thresholds
  cs_to_keep <-Reduce(intersect,list(good_lbf_cs, good_purity_mean_cs, good_purity_min_cs,good_pval_cs))

  susie_output_filtered <- susie_output

  # Print verbose information
  if(verbose){
    print(paste0(c("CS with lbf above threshold: ",good_lbf_cs),collapse = " "))
    print(paste0(c("CS with mean purity above threshold: ",good_purity_mean_cs),collapse = " "))
    print(paste0(c("CS with min purity threshold: ",good_purity_min_cs),collapse = " "))
    print(paste0(c("P-values of the index SNP for each CS: ",min_pval_per_cs),collapse = " "))
    print(paste0(c("CS with index SNP P-value below threshold: ",good_pval_cs),collapse = " "))
    print(paste0(c("CS which sutisfy all thresholds: ",cs_to_keep),collapse = " "))
  }

  # If there is at least one CS that passed threshold - report it, otherwise return NULL
  if(length(cs_to_keep)!=0){

    susie_output_filtered$sets$cs <- susie_output_filtered$sets$cs[paste0("L",cs_to_keep)]
    susie_output_filtered$sets$purity <- susie_output_filtered$sets$purity[paste0("L",cs_to_keep),]
    susie_output_filtered$sets$coverage <- susie_output_filtered$sets$coverage[which(susie_output_filtered$sets$cs_index %in% cs_to_keep)]
    susie_output_filtered$sets$cs_index <- cs_to_keep

    return(susie_output_filtered)

  } else {
    print("None of the credible sets satisfy selected thresholds")
    return(NULL)
  }

}
