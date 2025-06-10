#' Calculate Pairwise Colocalization and Adjacency Matrices
#'
#' @description
#' This function performs pairwise colocalization across signals in the same
#' clique and computes an adjacency matrix from approximated Bayes factors (ABFs).
#' It iterates through all unique pairs of cliques and applies the `icoloc`
#' function to determine the posterior probability of signals
#' sharing a common causal variant. An adjacency matrix is then
#' derived by thresholding these probabilities.
#'
#' @param clique.abfs A list of numeric vectors, where each vector contains
#'   approximated Bayes factors (ABFs) for SNPs within a specific clique.
#'   The names of the elements within each vector should correspond to the SNP identifiers.
#' @param threshold A numeric value between 0 and 1 (inclusive) representing the
#'   colocalization probability threshold. If the PP.H4.abf between two cliques
#'   is greater than or equal to this threshold, their corresponding entry in the
#'   adjacency matrix will be set to 1; otherwise, it will be 0. Defaults to 0.75.
#'
#' @return A list containing two matrices:
#'   \item{coloc.matrix}{A symmetric numeric matrix where `coloc.matrix[i, j]`
#'     represents the PP.H4.abf (posterior probability of shared causal variant)
#'     between `clique.abfs[[i]]` and `clique.abfs[[j]]`.}
#'   \item{adiacency.matrix}{A symmetric binary matrix (0s and 1s) representing
#'     the adjacency based on the `threshold`. `adiacency.matrix[i, j]` is 1 if
#'     the colocalization probability between `clique.abfs[[i]]` and `clique.abfs[[j]]`
#'     is at or above the `threshold`, and 0 otherwise. The diagonal elements are set to 0.}
#'
#' @details
#' The function uses `icoloc` for the colocalization analysis, which requires
#' input dataframes with 'snp' and 'lABF' columns. The parameters `p1`, `p2`, and `p12`
#' within `icoloc` are fixed at 1e-4, 1e-4, and 1e-5 respectively, representing
#' prior probabilities for different colocalization hypotheses.
#' The `adiacency.matrix` is derived from the `coloc.matrix` by applying the
#' specified `threshold`. Only the lower triangle of `coloc.matrix` is explicitly
#' calculated due to symmetry, and then mirrored.
#'
#' @examples
#' # Example usage
#' # For demonstration, create some dummy ABF lists
#' set.seed(123)
#' clique1_abfs <- setNames(rnorm(10, 5, 1), paste0("snp", 1:10))
#' clique2_abfs <- setNames(rnorm(10, 4, 1), paste0("snp", 6:15))
#' clique3_abfs <- setNames(rnorm(10, 6, 1), paste0("snp", 11:20))
#'
#' my_clique_abfs <- list(clique1_abfs, clique2_abfs, clique3_abfs)
#'
#' # Calculate colocalization and adjacency matrices
#' results <- coloc.abf.matrix(clique.abfs = my_clique_abfs, threshold = 0.8)
#'
#' # View the colocalization matrix
#' print(results$coloc.matrix)
#'
#' # View the adjacency matrix
#' print(results$adiacency.matrix)
#'
#' # Example with a lower threshold
#' results_low_thresh <- coloc.abf.matrix(my_clique_abfs, threshold = 0.5)
#' print(results_low_thresh$adiacency.matrix)
#'
#' @importFrom Matrix tril
#' @export
#'
coloc.abf.matrix <- function(clique.abfs=list(),threshold=0.75){
  W_comm=matrix(nrow = length(clique.abfs),ncol=length(clique.abfs))
  for(i in 1:(nrow(W_comm)-1)){
    for(j in (i+1):nrow(W_comm)){
      res=icoloc(
        df1=data.frame(snp=names(clique.abfs[[i]]),lABF=unlist(clique.abfs[[i]])),
        df2=data.frame(snp=names(clique.abfs[[j]]),lABF=unlist(clique.abfs[[j]])),
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = 1e-5
      )
      W_comm[i,j]<-W_comm[j,i]<-res$summary$PP.H4.abf
    }
  }
  adiacency.matrix=Matrix::tril(W_comm)
  diag(adiacency.matrix)=0
  adiacency.matrix[adiacency.matrix<threshold]=0
  adiacency.matrix[adiacency.matrix>0]=1

  res=list(coloc.matrix=W_comm,adiacency.matrix=adiacency.matrix)
  return(res)
}
