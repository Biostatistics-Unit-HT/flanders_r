#' Coloc percolation engine
#'
#' @description
#' This function takes fine-mapping data (approximated Bayes factors) for multiple
#' traits and performs a "colocalization-guided" percolation analysis.
#' It first calculates pairwise colocalization probabilities (PP.H4.abf) between
#' all input traits/regions, then thresholds these probabilities to create a binary
#' adjacency matrix. This matrix is then used to identify fully connected components
#' or "cliques" of colocalized signals. Optionally, it can apply `icoloc.regularization`
#' to the fine-mapping data before colocalization.
#'
#' @param finemap.data = matrix of log bayes factors obtained by finemapping tool (i.e. Susie), stored in assay(sce, "X"). Each column should correspond to a credible set while each row is a different SNP
#' @param pp.H4.threshold = numeric from 0 to 1 PP H4 threshold to call positive colocalization. Default is set to 0.75
#' @param icoloc.reg A logical value. If `TRUE`, the `icoloc.regularization`
#'   function will be applied to each column of `finemap.data` prior to
#'   calculating colocalization probabilities. Defaults to `TRUE`.
#'
#' @return
#' A list object with the following elements:
#'     \describe{
#'\item{ cp.object }{list object with the following items:
#'    \itemize{
#'      \item{indexes of each member of each community}
#'      \item{labels for the members in each community}
#'      \item{index of shared nodes}
#'      \item{index of shared nodes}
#'      \item{shared nodes labels}
#'      \item{colocalization matrix of final communities}
#'    }
#'}
#'.    \item{ abf.list }{list with vector of lABF per final community}
#'     \item{ coloc.matrix }{ adiacency matrix of colocalising credible sets (usefull for downstream plotting)}
#' }
#' @export
#'
#'@import dplyr qgraph coloc
#'
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

    names(results_unweighted) <- c("list.of.communities.numbers","list.of.communities.labels",
                                   "shared.nodes.numbers","shared.nodes.labels",
                                   "coloc.matrix")

    cl.list=list(list())
    cl.list[[1]]$clique.map=results_unweighted$list.of.communities.numbers

    fin.abf.list=extract.abf.list.matrix(data=finemap.data,clique.list =cl.list )

    returned_object <- c(results_unweighted,list( method = "coloc"))
    returned_object=list(cp.object=returned_object,abf.list=fin.abf.list,coloc.matrix=coloc.matrix)

    return(returned_object)

  } else {
    print("Only one clique detected")
    joint.finemapping=list(rowSums(finemap.data))
    returned_object=list(cp.object=NULL,abf.list=joint.finemapping,coloc.matrix=coloc.matrix)
    
    return(returned_object)
  }
}
