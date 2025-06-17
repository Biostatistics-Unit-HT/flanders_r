#' Perform Recursive Colocalization Percolation on a Matrix
#'
#' @description
#' This function takes an initial set of cliques (communities of colocalized
#' signals) and iteratively refines them through a recursive colocalization
#' percolation process. In each iteration, it aggregates the approximated Bayes factors
#' (ABFs) for the SNPs within the current level's communities, performs a
#' pairwise colocalization analysis between these aggregated communities, and then
#' identifies new, higher-level communities based on the colocalization results.
#' This process continues until the community structure stabilizes.
#' It also identifies shared nodes (signals) between the final communities.
#'
#' @param W An adjacency matrix (or a `qgraph` object from which `getWmat` can
#'   extract an adjacency matrix) representing initial connections between traits/regions.
#'   This is typically the thresholded colocalization matrix from a previous step.
#' @param cliques A list of numeric vectors, where each vector contains the
#'   indices of traits/regions (columns in `finemap.data`) belonging to an
#'   initial set of cliques. This is often the output from an initial unweighted
#'   clique detection on the first colocalization matrix.
#' @param labels A character vector containing the labels (names) corresponding
#'   to the traits/regions in `finemap.data`. These labels are used for the
#'   output communities.
#' @param finemap.data A matrix where rows represent SNPs and columns represent
#'   different traits or genomic regions. The values within the matrix should be
#'   approximated Bayes factors (lABF values) for each SNP in each trait/region.
#'
#' @return A list containing the results of the colocalization percolation:
#'   \item{list.of.communities.numbers}{A list of numeric vectors, where each
#'     vector contains the final numeric indices of traits/regions belonging
#'     to each detected community after percolation.}
#'   \item{list.of.communities.labels}{A list of character vectors, where each
#'     vector contains the labels (names) of traits/regions belonging to each
#'     detected community.}
#'   \item{shared.nodes.numbers}{A numeric vector of indices of traits/regions
#'     that are shared between two or more of the final communities. Empty if no
#'     shared nodes are found or if there's only one community.}
#'   \item{shared.nodes.labels}{A character vector of labels for traits/regions
#'     that are shared between two or more of the final communities. Empty if no
#'     shared nodes are found or if there's only one community.}
#'   \item{coloc.matrix}{The final binary adjacency matrix (`W_comm`) derived
#'     from the colocalization analysis of the aggregated ABFs at the last
#'     percolation level. This matrix reflects the connections between the
#'     detected communities.}
#'
#' @details
#' The function operates recursively until the community structure stabilizes:
#' \enumerate{
#'   \item **Initialize:** It starts with the provided `cliques` as the first
#'     level of communities.
#'   \item **Aggregate ABFs:** It calls `extract.abf.list.matrix` to combine
#'     the ABFs for SNPs within each community at the current level.
#'   \item **Colocalization Between Communities:** It then uses `coloc.abf.matrix`
#'     to perform pairwise colocalization analysis between these newly aggregated
#'     ABF sets. This results in a new adjacency matrix (`W_comm`) where nodes
#'     are now the communities from the previous level.
#'   \item **Identify New Communities:** `qgraph::qgraph` and `unweighted` (from `qgraph` or `igraph`)
#'     are used to detect fully connected components within this new adjacency matrix.
#'     These form the communities for the next level. `igraph::components`
#'     is used to check if the number of communities has stabilized.
#'   \item **Recursion:** Steps 2-4 are repeated. If the number of detected
#'     communities remains the same across iterations, the process stops.
#'   \item **Unravel Hierarchy:** Once stabilization is achieved, the function
#'     unravels the hierarchical `structure_level` to determine the final,
#'     most granular set of trait/region memberships for each community.
#'   \item **Identify Shared Nodes:** It then finds any traits/regions that are
#'     common to more than one of the final communities.
#'   \item **Assign Labels:** Finally, it translates the numeric community and
#'     shared node lists into human-readable labels using the input `labels` vector.
#' }
#' This percolation approach allows for the discovery of nested or overlapping
#' colocalization patterns, building from initial trait-level colocalizations
#' to higher-order community structures.
#'
#' @seealso
#' \code{\link[qgraph]{qgraph}}, \code{\link[qgraph]{getWmat}},
#' \code{\link[igraph]{graph_from_adjacency_matrix}}, \code{\link[igraph]{components}}
#'
#' @importFrom qgraph qgraph getWmat
#' @importFrom igraph graph_from_adjacency_matrix components
#' @keywords internal
#'
results_cp_coloc_matrix <- function(W=Wmat,
                                    cliques=cliques_unweighted,
                                    labels=labels,
                                    finemap.data=finemap.data
){

  #retrieve abfs
  structure_level=list()
  structure_level[[1]]=list(clique.map=NA)
  structure_level[[1]]$clique.map=cliques
  level=1

  tmp.abf=extract.abf.list.matrix(data=finemap.data,clique.list=structure_level)
#  print(6)

  W_comm=coloc.abf.matrix(clique.abfs = tmp.abf)
  W_comm=W_comm$adiacency.matrix
  fully.con=unweighted(qgraph::getWmat(qgraph::qgraph(W_comm)))
  W_i_comm <- igraph::graph_from_adjacency_matrix(as.matrix(W_comm),mode = "max")
  members <- igraph::components(W_i_comm)$membership
#  print(7)

  while(length(fully.con)!=length(unique(members))){
    print(level)
    level=level+1
    structure_level[[level]]=list(clique.map=NA)
    structure_level[[level]]$clique.map=fully.con
    tmp.abf=extract.abf.list.matrix(data=finemap.data, clique.list=structure_level)

    W_comm=coloc.abf.matrix(tmp.abf)
    W_comm=W_comm$adiacency.matrix
    fully.con=unweighted(qgraph::getWmat(qgraph::qgraph(W_comm)))
    W_i_comm <- igraph::graph_from_adjacency_matrix(as.matrix(W_comm))
    members <- (igraph::components(W_i_comm)$membership)
    if(length(fully.con)==length(unique(members))){
      level=level+1
      structure_level[[level]]=list(clique.map=NA)
      structure_level[[level]]$clique.map=fully.con

    }
  }

  final.groups=last(structure_level)$clique.map
  if(level>1){
    for( livello in length(structure_level):2){
      for( i in 1:length(final.groups)){

        final.groups[[i]]=unique(unlist(structure_level[[livello-1]]$clique.map[final.groups[[i]]]))

      }
    }
  } else {

    split=c()
    for(i in 1:length(fully.con)){
      split[[i]]=cliques[fully.con[[i]]]
    }

    #members <- igraph::components(W_i_comm)$membership
    #split <- split(cliques_unweighted, fully.con)
    final.groups <- lapply(split, function(x) sort(unique(unlist(x))))
  }

  message("adiacency created")
  communities=final.groups

  #if there is more than one community...
  #get shared nodes
  #if there are no shared nodes, create empty variable

  if (length(communities) > 1) {
    shared <- list()
    count <- 1
    for (i in 1:(length(communities) - 1)) {
      for (j in (i + 1):length(communities)) {
        shared[[count]] <- Reduce(intersect, list(communities[[i]],communities[[j]]))
        count <- count + 1
      }
    }
    shared <- sort(unique(unlist(shared)))
    if (length(shared) == 0) {shared <- c()}
  }

  #if there are communities...
  #create list of communities with labels from qgraph object instead of node numbers
  ### Add fix here

  if (length(communities) > 0) {
    communities_labels <- vector("list", length(communities))
    for (i in 1:length(communities)) {
      communities_labels[[i]] <- labels[communities[[i]]]
    }
  }

  #if there are no communities...
  if (length(communities) == 0) {
    communities_labels <- list()
  }

  #if there is more than one community...
  #create vector with shared nodes with labels from qgraph object instead of node numbers
  if (length(communities) > 1) {
    shared_labels <- labels[shared]
  }
  #if there is less than two communities
  if (length(communities) < 2) {
    shared <- c()
    shared_labels <- c()
  }
  #if there are no shared nodes, create empty variable
  if (length(shared_labels) == 0) {shared_labels <- c()}

  #return community lists, shared nodes vectors, and isolated nodes vectors
  return(list(list.of.communities.numbers=communities,
              list.of.communities.labels=communities_labels,
              shared.nodes.numbers=shared,
              shared.nodes.labels=shared_labels,
              coloc.matrix=W_comm))
}
