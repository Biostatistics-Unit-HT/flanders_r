# Sodbo Sharapov, Human Technopole (c) 2024

#' Meta-analysis of PIPs for multiple credible sets
#'
#' @description
#' This function inputs the list of vectors, each one containing PIPs for same
#' sets of SNPs and outputs a single vector with meta-analyzed PIPs.
#'
#' @param list_pips A SuSiE_R output (e.g. susieR::fitted_rss()).
#' @param verbose A boolean value. If TRUE, the log information will be printed.
#'
#' @return A single numveric vector of PIPs
#'
#' @details
#' This function inputs the list of vectors, each one containing PIPs for same
#' sets of SNPs and outputs a single vector with meta-analyzed PIPs.
#'
#' @examples
#' \dontrun{
#' pip.ma <- pip.meta(
#'     list_pips = susie_output,
#'     verbose = TRUE
#' )
#' }
#' @export
pip.meta <- function(
    list_pips = NULL,
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
