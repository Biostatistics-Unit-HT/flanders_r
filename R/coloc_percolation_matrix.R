#' Coloc Percolation Engine
#'
#' @description
#' This function takes fine-mapping data (approximated Bayes factors) for multiple
#' traits and performs a colocalization-guided percolation analysis.
#' It first calculates pairwise colocalization probabilities (PP.H4.abf) between
#' all input traits or regions, then thresholds these probabilities to create a binary
#' adjacency matrix. This matrix is used to identify fully connected components
#' (cliques) of colocalized signals. Optionally, it can apply `icoloc.regularization`
#' to the fine-mapping data before colocalization.
#'
#' @param finemap.data A matrix of log Bayes factors obtained from a fine-mapping tool
#'   (e.g., SuSiE), typically stored in `assay(sce, "X")`. Each column corresponds to a
#'   credible set and each row to a SNP.
#' @param pp.H4.threshold Numeric between 0 and 1. Posterior probability threshold (PP.H4)
#'   to call colocalization. Default is 0.75.
#' @param icoloc.reg Logical. If `TRUE`, the `icoloc.regularization` function is applied
#'   to each column of `finemap.data` before colocalization. Default is `TRUE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{cp.object}{List containing:
#'     \itemize{
#'       \item{Indexes of members in each community}
#'       \item{Labels of community members}
#'       \item{Indexes of shared nodes}
#'       \item{Labels of shared nodes}
#'       \item{Colocalization matrix for final communities}
#'     }
#'   }
#'   \item{abf.list}{List of vectors with log Bayes factors per final community}
#'   \item{coloc.matrix}{Adjacency matrix of colocalizing credible sets (useful for downstream plotting)}
#'   \item{post_cp_tmp_module}{Data frame with `cs_name` and a temporary module name (random string)}
#' }
#'
#' @import dplyr qgraph coloc
#' @export
#'
coloc_percolation_matrix <- function(finemap.data = cs.data,
                                     pp.H4.threshold=0.75,
                                     icoloc.reg=TRUE){

  ### Step 1 create coloc matrix
  if(icoloc.reg==TRUE){
    finemap.data=apply(finemap.data,2,icoloc.regularization)
  }
  coloc.matrix=matrix(nrow=ncol(finemap.data),ncol=ncol(finemap.data))
  rownames(coloc.matrix)<-colnames(coloc.matrix)<-colnames(finemap.data)
  col.name=colnames(coloc.matrix)
  pairs=expand.grid(col.name,col.name,stringsAsFactors = F)
  pairs <- unique(t(apply(pairs, 1, sort)))
  keep=apply(pairs,1,function(x)x[1]!=x[2])
  pairs=pairs[keep,]

  for(i in 1:nrow(pairs)){
    res=icoloc.ht(df1=data.frame(snp=row.names(finemap.data),lABF=finemap.data[,pairs[i,1]]),
                  df2=data.frame(snp=row.names(finemap.data),lABF=finemap.data[,pairs[i,2]]))

    coloc.matrix[pairs[i,1],pairs[i,2]]=res$summary$PP.H4.abf
  }

  ## prune all values below colocalization threshold
  coloc.matrix[coloc.matrix<pp.H4.threshold]=0
  coloc.matrix[coloc.matrix>0]=1

  # create colocalization graph
  W=qgraph::qgraph(coloc.matrix,labels=colnames(coloc.matrix),DoNotPlot=T)
  trait.list=colnames(coloc.matrix)
  Wmat <- qgraph::getWmat(W)
  labels <- as.vector(W$graphAttributes$Nodes$labels)

  # identify fully connected cliques
  cliques_unweighted <- unweighted(Wmat)

  if(length(cliques_unweighted)>1){

    ## perform recursive coloc percolation
    results_unweighted <- results_cp_coloc_matrix(W = Wmat,
                                                  cliques = cliques_unweighted,
                                                  labels =  labels,
                                                  finemap.data =finemap.data
    )

    cl.list=list(list())
    cl.list[[1]]$clique.map=results_unweighted$list.of.communities.numbers

    fin.abf.list=extract.abf.list.matrix(data=finemap.data,clique.list =cl.list )

    returned_object <- c(results_unweighted,list( method = "coloc"))
    returned_object=list(
      cp.object=returned_object,
      abf.list=fin.abf.list,
      coloc.matrix=coloc.matrix)

    # Add to object also list of cs with random string module label
    returned_object$post_cp_tmp_module <- assign_post_cp_tmp_module_names(returned_object, trait.list)
    return(returned_object)

  } else {
    print("Only one clique detected")
    joint.finemapping=list(rowSums(finemap.data))
    returned_object=list(
      cp.object=NULL,
      abf.list=joint.finemapping,
      coloc.matrix=coloc.matrix)
    # Add to object also list of cs with random string module label
    returned_object$post_cp_tmp_module <- assign_post_cp_tmp_module_names(returned_object, trait.list)
    return(returned_object)
  }
}
