#' Extract and Process lABF Matrix from Annotated Data Across Chromosomes
#'
#' @description
#' This function extracts a matrix of log-approximated Bayes factors (lABFs)
#' from a `SingleCellExperiment` (SCE) object, specifically for a set of
#' credible sets (CSs) provided in `cs.list`. It supports multiple chromosomes
#' by processing each chromosome separately and then combining the results.
#' For each CS, any zero lABF values are replaced by the column-specific
#' `min_res_labf` to ensure imputed background signal.
#'
#' @param sce.obj A `SingleCellExperiment` object containing SNP annotations and assay data.
#' @param cs.list A character vector of credible set names (column names in `sce.obj`) to extract.
#'
#' @return A matrix of lABF values, where columns represent credible sets and rows represent SNPs. 
#'         SNPs across multiple chromosomes are supported, and the matrix is padded with `NA` 
#'         where certain credible sets do not appear on a given chromosome.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **Subset by Credible Sets:** Filters the `sce.obj` to include only the credible sets specified in `cs.list`.
#'   \item **Loop Over Chromosomes:** Iterates over each chromosome present in the subset.
#'   \item **Define Genomic Boundaries per Chromosome:** For each chromosome, defines boundaries using the min `start` and max `end` positions.
#'   \item **Filter SNPs Within Bounds:** Filters rows (SNPs) to retain those within the chromosome-specific genomic region.
#'   \item **Extract and Modify Matrix:** Extracts the lABF matrix from the 'X' assay. For each CS, replaces `0` values with `min_res_labf`.
#'   \item **Column Alignment Across Chromosomes:** Ensures all matrices have the same set of columns by padding with `NA` where needed.
#'   \item **Combine Results:** Row-binds the per-chromosome matrices into a single output matrix.
#' }
#'
#' @seealso
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @import zellkonverter SingleCellExperiment
#'
#' @keywords internal
#'
extract.ann.data.labf.matrix <- function(sce.obj = ann.data,
                                         cs.list = member.list) {
  sce.obj = sce.obj[, cs.list]
  
  chromosomes = unique(as.character(sce.obj$chr))
  all_matrices = list()
  
  for (chr in chromosomes) {
    chr_sce = sce.obj[, sce.obj$chr == chr]
    
    # Get chromosome-specific boundaries
    boundries = c(min(chr_sce$start), max(chr_sce$end))
    
    # Filter rows for the current chromosome and region
    idx = which(rowData(chr_sce)$chr == chr &
                  rowData(chr_sce)$pos >= boundries[1] &
                  rowData(chr_sce)$pos <= boundries[2])
    
    chr_sce = chr_sce[idx, ]
    
    if (nrow(chr_sce) == 0 || ncol(chr_sce) == 0) next  # Skip empty cases
    
    cs.matrix = chr_sce@assays@data$X
    colnames(cs.matrix) = colnames(chr_sce)
    rownames(cs.matrix) = rownames(chr_sce)
    cs.matrix = as.matrix(cs.matrix)
    
    for (i in colnames(cs.matrix)) {
      cs.matrix[which(cs.matrix[, i] == 0), i] = chr_sce[, i]$min_res_labf
    }
    
    all_matrices[[chr]] = cs.matrix
  }
  
  # Get all unique column names across matrices
  all_cols = unique(unlist(lapply(all_matrices, colnames)))
  
  all_matrices_padded = lapply(all_matrices, pad_matrix, all_cols = all_cols)
  final_matrix = do.call(rbind, all_matrices_padded)
  return(final_matrix)
}