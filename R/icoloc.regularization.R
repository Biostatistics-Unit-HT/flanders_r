#' Regularize Log-Approximated Bayes Factors (lABFs) for Colocalization
#'
#' @description
#' This function regularizes a vector of log-approximated Bayes factors (lABFs)
#' by identifying a 99% credible set based on posterior inclusion probabilities (PIPs).
#' SNPs falling outside this 99% credible set have their lABF values imputed
#' to the minimum lABF observed for all SNPs outside the credible set. This
#' regularization can help in downstream colocalization analyses by providing
#' a consistent background lABF for non-credible SNPs.
#'
#' @param cs.vector A numeric vector of log-approximated Bayes factors (lABFs)
#'   for SNPs within a single credible set or genomic region of a particular trait/study.
#'
#' @return A numeric vector of regularized lABFs, with the same length and order
#'   as the input `cs.vector`. SNPs within the 99% credible set retain their
#'   original lABFs, while those outside have their lABFs imputed.
#'
#' @details
#' The regularization process involves several steps:
#' \enumerate{
#'   \item **Calculate Posterior Inclusion Probabilities (PIPs):** The input
#'     `cs.vector` (lABFs) is converted into posterior probabilities using the
#'     `logbf_to_pp` function. These probabilities can be interpreted as PIPs
#'     for each SNP.
#'   \item **Identify 99% Credible Set:** SNPs are sorted by their PIPs in
#'     decreasing order. A cumulative sum of PIPs is then calculated, and the
#'     99% credible set is defined as the minimum number of top SNPs whose
#'     cumulative PIP exceeds 0.99.
#'   \item **Mark Credible SNPs:** A logical flag is assigned to each SNP,
#'     indicating whether it belongs to the 99% credible set (`TRUE`) or not (`FALSE`).
#'   \item **Impute lABFs:** For SNPs *outside* the 99% credible set (where the
#'     flag is `FALSE`), their lABF values are replaced with the minimum lABF
#'     value observed among all SNPs outside this 99% credible set. This imputation
#'     ensures that these SNPs contribute minimally but consistently to downstream
#'     analyses without introducing spurious signals from very small lABFs close to zero.
#'   \item **Restore Original Order:** Finally, the regularized lABFs are
#'     reordered back to their original input order based on their initial indices.
#' }
#' This regularization is often applied as a preprocessing step before performing
#' colocalization analysis with functions like `icoloc`, as it helps to focus
#' the analysis on the most relevant SNPs within each credible set.
#'
#' @seealso
#' \code{\link{logbf_to_pp}} (for converting lABFs to posterior probabilities)
#'
#' @keywords internal
#'
icoloc.regularization <- function(cs.vector){
  cs.pip=logbf_to_pp(cs.vector)
  cs.vector=cbind(1:length(cs.vector),cs.vector,cs.pip,FALSE)
  cs.vector=cs.vector[order(cs.vector[,3],decreasing=T),]
  last.idx=min(which(cumsum(cs.vector[,3])>0.99))
  cs.vector[1:last.idx,4]=1
  cs.vector[cs.vector[,4]==0,2]=min(cs.vector[cs.vector[,4]==0,2],na.rm=T)
  cs.vector=cs.vector[order(cs.vector[,1],decreasing=F),]
  cs.vector[,2]
}
