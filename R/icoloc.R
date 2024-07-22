#' Perform Colocalization Analysis Using lABF Values
#'
#' This function performs colocalization analysis between two datasets using lABF values.
#' It computes the posterior probability that each SNP is the causal variant for a shared signal.
#'
#' @param df1 A data frame for the first dataset containing SNPs and their lABF values.
#' @param df2 A data frame for the second dataset containing SNPs and their lABF values.
#' @param p1 A numeric value representing the prior probability for the first dataset. Default is 1e-4.
#' @param p2 A numeric value representing the prior probability for the second dataset. Default is 1e-4.
#' @param p12 A numeric value representing the prior probability for both datasets. Default is 1e-5.
#'
#' @return A list containing the colocalization summary and results by SNP.
#' @export
#'
#' @examples
#' \dontrun{
#' dataset1 <- data.frame(snp = c("rs1", "rs2", "rs3"), lABF = c(1.2, 2.3, 1.5))
#' dataset2 <- data.frame(snp = c("rs1", "rs2", "rs3"), lABF = c(1.5, 2.1, 1.7))
#' result <- icoloc.ht(dataset1, dataset2)
#' print(result$summary)
#' print(result$results)
#' }
icoloc.ht <- function(
    df1 = dataset1,
    df2 = dataset2,
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 1e-5
){

  df1 <- df1 %>% rename("lABF.df1" = "lABF")
  df2 <- df2 %>% rename("lABF.df2" = "lABF")

  p1 <- coloc:::adjust_prior(p1, nrow(df1), "1")
  p2 <- coloc:::adjust_prior(p2, nrow(df2), "2")

  merged.df <- merge(df1, df2, by = "snp")
  p12 <- coloc:::adjust_prior(p12, nrow(merged.df), "12")

  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)

  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps = common.snps, pp.abf)

  colo.res <- list(summary = results, results = merged.df, priors = c(p1 = p1, p2 = p2, p12 = p12))
  class(colo.res) <- c("coloc_abf", class(colo.res))

  ## Save coloc summary
  colo.sum <- data.frame(t(colo.res$summary))

  ## Save coloc result by SNP
  colo.full_res <- colo.res$results %>% dplyr::select(snp, lABF.df1, lABF.df2, SNP.PP.H4)

  ## Organise all in a list (composed of summary + results)
  coloc.final <- list(summary = colo.sum, results = colo.full_res)
  return(coloc.final)
}
