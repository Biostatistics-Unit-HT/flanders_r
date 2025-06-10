#' extract.ann.data.labf.matrix
#'
#' @param sce.obj = SingleCellExperiment object containing the CS information
#' @param cs.list = vector with names of the credible sets to be extracted
#'
#' @return matrix with lABFs values. Columns are credible sets IDs while rows are SNPs.
#'
#' @export
#'
#'@import zellkonverter SingleCellExperiment
#'
#'
extract.ann.data.labf.matrix=function(sce.obj=ann.data
                                      ,cs.list=member.list ){

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
