#' Identify Simple Colocalization Modules and Assess Connectivity
#'
#' @description
#' This function takes a data frame of pairwise colocalization results, filters
#' them based on a posterior probability of shared causal variant (PP.H4.abf)
#' threshold, and then defines modules (communities) of colocalizing genomic
#' regions or traits using a graph-based approach. It also assesses whether
#' these identified modules are fully connected, categorizing them into
#' partly connected and fully connected.
#'
#' @param coloc.results A data frame containing pairwise colocalization results.
#'   It is expected to have at least the following columns:
#'   \itemize{
#'     \item `PP.H4.abf`: Numeric, the posterior probability of shared causal variant.
#'     \item `hit1`: Character, identifier for the first genomic region/trait in the pair.
#'     \item `hit2`: Character, identifier for the second genomic region/trait in the pair.
#'   }
#' @param pp.H4.threshold A numeric value (between 0 and 1) representing the
#'   threshold for `PP.H4.abf`. Only pairs with a `PP.H4.abf` strictly greater
#'   than this threshold will be considered for module definition. Defaults to 0.75.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item `modules`: A `data.table` with two columns:
#'       \itemize{
#'         \item `cs`: Character, representing the names of the genomic regions/traits.
#'         \item `module`: Numeric, the assigned module ID for each `cs`.
#'       }
#'       The data.table is ordered by module ID.
#'     \item `partly_connected_modules`: A numeric vector of module IDs that are
#'       not fully connected (i.e., not all traits within the module colocalize
#'       with every other trait in that same module).
#'     \item `fully_connected_modules`: A numeric vector of module IDs that are
#'       fully connected.
#'   }
#'   If no pair-wise colocalization results meet the `pp.H4.threshold`, the
#'   function returns `NULL` and prints a message.
#'
#' @details
#' The function first filters the input `coloc.results` to retain only those
#' pairs where `PP.H4.abf` exceeds `pp.H4.threshold`.
#'
#' An `igraph` graph is then constructed where nodes are the unique `hit1` and `hit2`
#' identifiers (genomic regions/traits) and edges exist between colocalizing pairs.
#' Connected components within this graph are identified as `modules`.
#'
#' For each identified module, the function checks if all nodes within that
#' module are fully connected (i.e., if every trait within the module colocalizes
#' with every other trait in the same module). It calculates the `expected_edges`
#' for a fully connected subgraph of `n` nodes as `n * (n - 1) / 2`.
#' Modules that are not fully connected might indicate the need for further
#' hierarchical colocalization analysis or suggest complex relationships.
#'
#' Progress messages are printed periodically during module assessment.
#'
#' @importFrom igraph graph_from_data_frame components V subgraph ecount
#'
#' @export
#'
simple_modules_connectivity <- function(coloc.results, pp.H4.threshold=0.75) {  
  
  ## Filter pair-wise coloc results according to PPH4 threshold
  coloc.results <- coloc.results[round(as.numeric(coloc.results$PP.H4.abf), 2) > pp.H4.threshold, ]

  if(nrow(coloc.results)>0){

    ## Define modules based of colocalisation results
    coloc.graph <- igraph::graph_from_data_frame(coloc.results[, c("hit1","hit2")])
    modules <- igraph::components(coloc.graph)
    modules <- data.frame(cs = names(modules$membership), module = modules$membership)
    modules <- modules[order(modules$module), ]
    #modules.size <- table(modules$module)
    
    ## Identify modules that need further refining (not fully connected modules!)
    partly_connected <- c()
    fully_connected <- c()
    
    for (n_module in seq_along(unique(modules$module))) {
      traits <- modules$cs[modules$module == n_module]
      a <- igraph::V(coloc.graph)
      tmp <- subgraph(coloc.graph, which(names(a) %in% traits))
      
      expected_edges <- (length(traits)^2 - length(traits)) / 2
      
      # Check whether the subgraph is fully connected
      if (ecount(tmp) != expected_edges) {
        partly_connected <- c(partly_connected, n_module)
      } else {
        fully_connected <- c(fully_connected, n_module)
      }
      
      # Print every 100 iterations
      if (n_module %% 500 == 0 || n_module == length(unique(modules$module))) {
        message(sprintf("Processed %d of %d modules...", n_module, length(unique(modules$module))))
      }
    }
  return(
    list(partly_connected_modules = partly_connected,
      fully_connected_modules = fully_connected)
  )

  } else {
    message("No pairwise colocalization test found with PPH4 >", pp.H4.threshold)
    return(NULL)
  }
}
