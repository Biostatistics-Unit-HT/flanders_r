#' Extract and Process lABF Matrix from Annotated Data
#'
#' @description
#' This function extracts a matrix of log-approximated Bayes factors (lABFs)
#' from a `SingleCellExperiment` (SCE) object, specifically for a set of
#' credible sets (CSs) provided in `cs.list`. It filters the SCE object
#' to include only SNPs within the genomic boundaries covered by these CSs
#' and then imputes any zero lABF values for non-colocalizing SNPs with a
#' predefined minimum residual lABF from the SCE object's metadata.
#'
#' @param sce.obj = SingleCellExperiment object containing the CS information
#' @param cs.list = vector with names of the credible sets to be extracted
#'
#' @return matrix with lABFs values. Columns are credible sets IDs while rows are SNPs.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **Subset by Credible Sets:** The `sce.obj` is initially subsetted to
#'     include only the columns (credible sets) specified in `cs.list`.
#'   \item **Define Genomic Boundaries:** It then determines the minimum 'start'
#'     and maximum 'end' genomic positions from the subsetted `sce.obj` to define
#'     a comprehensive genomic region of interest. The chromosome is also identified.
#'   \item **Filter SNPs:** SNPs (rows) in the `sce.obj` are filtered to retain
#'     only those that fall within the determined chromosome and genomic boundaries.
#'   \item **Extract lABF Matrix:** The lABF values are extracted from the 'X'
#'     assay of the filtered `sce.obj` and converted into a standard matrix,
#'     ensuring correct row and column names (SNP IDs and credible set names).
#'   \item **Impute Zero lABFs:** For each column (credible set) in the
#'     extracted lABF matrix, any lABF values that are exactly `0` are replaced
#'     with the corresponding `min_res_labf` value from the `sce.obj`'s column
#'     metadata for that specific credible set. This imputation step is crucial
#'     for downstream colocalization analyses, as it provides a background
#'     lABF for SNPs that do not show evidence of association.
#' }
#'
#' @seealso
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @import zellkonverter SingleCellExperiment
#'
#' @keywords internal
#'
extract.ann.data.labf.matrix <- function(sce.obj=ann.data,
                                         cs.list=member.list){

  sce.obj=sce.obj[,cs.list]
  boundries=c(min(sce.obj$start),max(sce.obj$end))
  chr=unique(sce.obj$chr)
  idx=which(rowData(sce.obj)$chr==chr
            & rowData(sce.obj)$pos>=boundries[1]
            & rowData(sce.obj)$pos<=boundries[2])
  sce.obj=sce.obj[idx,]
  cs.matrix=sce.obj@assays@data$X
  colnames(cs.matrix)=colnames(sce.obj)
  rownames(cs.matrix)=rownames(sce.obj)
  cs.matrix=as.matrix(cs.matrix)
  for(i in colnames(cs.matrix)){

    cs.matrix[which(cs.matrix[,i]==0),i]=sce.obj[,i]$min_res_labf

  }
  return(cs.matrix)
}
