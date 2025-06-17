#' Meta-analysis per trait using colocalization results
#'
#' For each pair of colocalized credible sets (with high posterior probability and from different studies, same trait),
#' performs meta-analysis using \code{cs_metanalysis}.
#'
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing summary statistics.
#' @param coloc_res A data.frame or data.table with colocalization results, including columns \code{t1}, \code{t2}, \code{t1_study_id}, \code{t2_study_id}, \code{PP.H4.abf}, \code{hit1}, and \code{hit2}.
#'
#' @return A list of meta-analysis results, one for each pair of colocalized credible sets.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' library(zellkonverter)
#' library(dplyr)
#'
#' ann.data <- readH5AD(
#'   file = "/ssu/bsssu/anndata_finemapping_repository/cardinal_UKBB_europeans_ELGH_pakistani_cis_finemapping_full_FLIP_FIXED_20250611.h5ad",
#'   verbose = TRUE,
#'   reader = "R"
#' )
#'
#' coloc_res <- fread("/ssu/bsssu/cardinal_finemap_coloc/anndata/anndata_020525_coloc/results/cardinal_UKBB_europeans_ELGH_pakistani_colocalization_table_all_merged.txt")
#' results <- meta_analysis_per_trait(ann.data, coloc_res)
#' }
#' @export
meta_analysis_per_trait <- function(ad_or_sce, coloc_res) {

    coloc_res_to_MA <- coloc_res %>%
        filter(t1 == t2 & t1_study_id != t2_study_id & PP.H4.abf >= 0.75) %>%
        select(hit1, hit2)

    print(paste0("There are ", nrow(coloc_res_to_MA), " pairs of colocalized credible sets of same traits across studies."))

    # Convert each row to a list of two elements: hit1 and hit2
    cs_lists <- apply(coloc_res_to_MA, 1, function(row) as.list(row))

    i <- 0

    results <- lapply(
        cs_lists,
        function(cs_pair) {
            res <- cs_metanalysis(
                ad_or_sce,
                cs.list = unlist(cs_pair)
            )
            i <<- i + 1
            print(paste0("Meta-analysis ", i, " of ", length(cs_lists), " completed."))
            return(res)
        }
    )

    return(results)
    
}