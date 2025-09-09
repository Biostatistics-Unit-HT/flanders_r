#' Expand SuSiE Credible Sets by Thresholding
#'
#' Expands each credible set (CS) from a fitted \code{susie} or \code{susie_rss} 
#' object by applying a thresholding rule to log Bayes factors.  
#' The expanded set is defined as the union of the original CS and 
#' SNPs identified via \code{\link{find_threshold}}.
#'
#' @param fitted A \code{susie} or \code{susie_rss} object, typically the 
#'   output of \code{susieR::susie_rss()}.
#'
#' @return A named list of integer vectors, one per credible set.  
#'   Each vector contains SNP indices for the union of the original and 
#'   expanded credible sets, with SNP IDs as names.
#'
#' @details 
#' The function retrieves original CSs using \code{susie_get_cs()}, 
#' then applies \code{find_threshold()} to the corresponding lBF values.  
#' The expanded CS is the union of both sets.
#'
#' @seealso \code{\link{find_threshold}}, \code{\link[susieR]{susie_get_cs}}
#' @export
#'
#' @examples
#' \dontrun{
#' fitted <- susieR::susie_rss(z, R, n = N)
#' expanded <- expand_cs(fitted)
#' names(expanded)
#' }
expand_cs <- function(fitted) {
    cs_list  <- fitted$sets$cs
    snp_names <- colnames(fitted$lbf_variable)
    
    expanded <- lapply(seq_along(cs_list), function(i) {
        idx_original <- cs_list[[i]]
        names(idx_original) <- snp_names[idx_original]  # name with SNP IDs
    
        lbf_vec <- fitted$lbf_variable[i, ] ### are sure that the index of credible sets is the same here?
        idx_expanded <- find_threshold(lbf_vec, original_idx = idx_original)
    
        # Union of both
        idx_union <- union(idx_original, idx_expanded)
        names(idx_union) <- snp_names[idx_union]
        
        idx_union
    })
    names(expanded) <- paste0("L", seq_along(cs_list))
    return(expanded)
}