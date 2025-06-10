#' Identify Maximal Cliques in an Unweighted Adjacency Matrix
#'
#' @description
#' This function takes an unweighted adjacency matrix and identifies all
#' maximal cliques within the graph it represents. A maximal clique is a
#' complete subgraph that cannot be extended by adding any more vertices
#' from the graph. The function uses the `igraph` package for graph
#' representation and clique detection.
#'
#' @param W_unweighted A symmetric numeric matrix representing an unweighted
#'   adjacency matrix. Typically, values are 0 (no edge) or 1 (edge exists).
#'   Non-binary numeric values will be treated as existence of an edge.
#'
#' @return A list of numeric vectors. Each vector in the list represents a
#'   maximal clique, containing the numeric indices of the nodes (vertices)
#'   that form that clique.
#'
#' @details
#' The process involves two main steps:
#' \enumerate{
#'   \item **Graph Conversion:** The input `W_unweighted` matrix is converted
#'     into an `igraph` object using `igraph::graph_from_adjacency_matrix`.
#'     The `mode` is set to "undirected" as adjacency matrices in colocalization
#'     contexts are typically symmetric, and `weighted = NULL` ensures it's
#'     treated as an unweighted graph for clique finding.
#'   \item **Maximal Clique Extraction:** The `igraph::max_cliques` function is
#'     then applied to the `igraph` object to find all maximal cliques. These are
#'     returned as a list, where each element is an `igraph.vs` object.
#'     `lapply(as.vector)` converts these `igraph.vs` objects into simple
#'     numeric vectors of vertex indices.
#' }
#' This function is crucial for identifying tightly connected groups of traits
#' or genomic regions that show strong evidence of colocalization.
#'
#' @seealso
#' \code{\link[igraph]{graph_from_adjacency_matrix}}, \code{\link[igraph]{max_cliques}},
#' \code{\link[igraph]{cliques}} (for finding cliques of a specific size `k`)
#'
#' @importFrom igraph graph_from_adjacency_matrix max_cliques
#'
#' @export
#'
unweighted <- function(W_unweighted){

  # transform to igraph object to use specific functions
  W_i_unweighted <- igraph::graph_from_adjacency_matrix(
    W_unweighted,
    mode = "undirected", weighted = NULL)

  # extract cliques
  #cliques_unweighted <- igraph::cliques(W_i_unweighted, min = k, max = k) %>% lapply(as.vector)
  cliques_unweighted <- igraph::max_cliques(W_i_unweighted)%>% lapply(as.vector)
  # return cliques
  return(cliques_unweighted)
}
