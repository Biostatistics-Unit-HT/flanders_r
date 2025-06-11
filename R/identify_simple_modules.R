#' Identify Simple Colocalization Modules and Assess Connectivity
#'
#' @description
#' This function takes a data frame of pairwise colocalization results, filters
#' them based on a posterior probability of shared causal variant (`PP.H4.abf`)
#' threshold, and then defines modules (communities) of colocalizing genomic
#' regions or traits using a graph-based approach. Optionally, it checks whether
#' each module is fully connected.
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
#' @param check_connectivity Logical. If `TRUE` (default), each identified module
#'   is checked for full connectivity. If `FALSE`, modules are returned without
#'   further connectivity checks.
#'
#' @return If `check_connectivity = FALSE`, returns a `data.table` with:
#'   \itemize{
#'     \item `cs_name`: Character, names of genomic regions or traits.
#'     \item `module`: Integer, assigned module ID.
#'   }
#'   If `check_connectivity = TRUE`, returns a named list with two `data.tables`:
#'   \itemize{
#'     \item `partly_connected_modules`: Subset of modules that are not fully connected.
#'     \item `fully_connected_modules`: Subset of modules that are fully connected.
#'   }
#'   If no results pass the threshold, returns `NULL`.
#'
#' @details
#' The function first filters the input `coloc.results` to retain only those
#' pairs where `PP.H4.abf` exceeds the specified threshold.
#'
#' An `igraph` graph is built using the `hit1` and `hit2` columns as edges.
#' Connected components in this graph define modules. If connectivity is checked,
#' each module is evaluated to determine if all traits within it form a fully
#' connected clique (i.e., every node has an edge to every other node).
#'
#' Progress updates are printed every 500 modules during connectivity checks.
#'
#' @importFrom igraph graph_from_data_frame components V subgraph ecount
#' @import data.table
#'
#' @export
#'
identify_simple_modules <- function(coloc.results, pp.H4.threshold=0.75, check_connectivity=TRUE) {  
  
  ## Filter pair-wise coloc results according to PPH4 threshold
  coloc.results <- coloc.results[round(as.numeric(coloc.results$PP.H4.abf), 2) > pp.H4.threshold, ]

  if(nrow(coloc.results)>0){

    ## Define modules based of colocalisation results
    coloc.graph <- igraph::graph_from_data_frame(coloc.results[, c("hit1","hit2")])
    modules <- igraph::components(coloc.graph)
    modules <- data.table(cs_name = names(modules$membership), module = modules$membership)
    modules <- modules[order(modules$module), ]
    #modules.size <- table(modules$module)
    message(length(unique(modules$module)), " simple modules identified")
    
    if(check_connectivity==FALSE){ # Return modules as they are
      return(modules)
    }
    
    ## Identify modules that need further refining (not fully connected modules!)
    partly_connected <- c()
    fully_connected <- c()
    
    for(n_module in seq_along(unique(modules$module))) {
      traits <- modules$cs_name[modules$module == n_module]
      a <- igraph::V(coloc.graph)
      tmp <- subgraph(coloc.graph, which(names(a) %in% traits))
      
      expected_edges <- (length(traits)^2 - length(traits)) / 2
      
      # Check whether the subgraph is fully connected
      if (ecount(tmp) != expected_edges) {
        partly_connected <- c(partly_connected, n_module)
      } else {
        fully_connected <- c(fully_connected, n_module)
      }
      
      # Print every 500 iterations
      if (n_module %% 500 == 0 || n_module == length(unique(modules$module))) {
        message(sprintf("Processed %d of %d modules...", n_module, length(unique(modules$module))))
      }
    }
    partly_connected_df <- modules[module %in% partly_connected]
    fully_connected_df <- modules[module %in% fully_connected]
    
    return(
      list(partly_connected_modules = partly_connected_df,
        fully_connected_modules = fully_connected_df)
    )

  } else {
    message("No pairwise colocalization test found with PPH4 >", pp.H4.threshold)
    return(NULL)
  }
}
