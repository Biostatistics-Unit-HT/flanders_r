#' Assign Random Module Names to Post-colocalization Modules
#'
#' This function generates unique random identifiers for each post-colocalization module
#' found in a `coloc_percolation_matrix` output object. It processes each element in 
#' `cp.object$list.of.communities.labels` object and assigns a 30-character alphanumeric
#' module name. If no communities are present, it creates a single module from the `traits` object.
#'
#' @param groups A list object that contains a component `cp.object$list.of.communities.labels`,
#'   which should be a list of community labels (e.g., vectors or data frames).
#' @param traits A list of credible set names
#'
#' @return A `data.frame` with at least two columns:
#' \describe{
#'   \item{cs_name}{Original community structure or trait information.}
#'   \item{module}{Random 30-character alphanumeric module identifier.}
#' }
#'
#' @importFrom stringi stri_rand_strings
#' @keywords internal
#'
assign_post_cp_tmp_module_names <- function(groups, traits){
  all_modules <- list()
  if(length(groups$cp.object$list.of.communities.labels)>0){
    for(j in 1:length(groups$cp.object$list.of.communities.labels)){
      tmp <- data.frame(
        cs_name = groups$cp.object$list.of.communities.labels[[j]],
        module = stringi::stri_rand_strings(n=1, length=30, pattern="[A-Za-z0-9]")
      )
      all_modules <- rbind(all_modules, tmp)
    }
    
  } else {
    traits$module <- stringi::stri_rand_strings(n=1, length=30, pattern="[A-Za-z0-9]")
    all_modules <- rbind(all_modules, traits)
  }
  return(all_modules)
}
