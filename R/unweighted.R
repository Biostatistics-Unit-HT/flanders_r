




unweighted <- function(W_unweighted){

  #transform to igraph object to use specific functions
  W_i_unweighted <- igraph::graph_from_adjacency_matrix(
    W_unweighted,
    mode = "undirected", weighted = NULL)

  #extract cliques
  #cliques_unweighted <- igraph::cliques(W_i_unweighted, min = k, max = k) %>% lapply(as.vector)
  cliques_unweighted= igraph::max_cliques(W_i_unweighted)%>% lapply(as.vector)
  #return cliques
  return(cliques_unweighted)
}
