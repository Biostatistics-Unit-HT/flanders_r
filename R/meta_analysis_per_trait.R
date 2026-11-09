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
#' @import igraph
#' @import dplyr
#' @export
meta_analysis_per_trait <- function(ad_or_sce, coloc_res) {

    coloc_res_to_MA <- coloc_res %>%
        filter(t1 == t2 & t1_study_id != t2_study_id & PP.H4.abf >= 0.75)

    print(paste0("There are ", nrow(coloc_res_to_MA), " pairs of colocalized credible sets of same traits across studies."))

    modules <- identify_simple_modules(
        coloc.results = coloc_res_to_MA,
        check_connectivity = TRUE
        )

    print(paste0("There are ", length(unique(modules$partly_connected_modules$module)), " partly connected modules and ", length(unique(modules$fully_connected_modules$module)), " fully connected modules."))
    print(paste0("So ", length(unique(modules$partly_connected_modules$module)), " traits require coloc perculation and ", length(unique(modules$fully_connected_modules$module)), " traits can be meta-analyzed directly."))

    master.table.meta <- colData(ann.data) %>% as.data.table

    all_cs <- c(modules$partly_connected_modules$cs_name,modules$fully_connected_modules$cs_name)
    
    # Now master.table.meta should only contain rows with cs_name which will participate in the meta-analysis
    master.table.meta <- master.table.meta[master.table.meta$cs_name %in% all_cs,]

    all_mod_tab <- Reduce(rbind, modules)

    master.table.meta$module <- all_mod_tab$module[match(master.table.meta$cs_name,all_mod_tab$cs_name)]
    master.table.meta$module <- all_mod_tab$module[match(master.table.meta$cs_name,all_mod_tab$cs_name)]
    master.table.meta <- master.table.meta[order(master.table.meta$module),]

    modules_for_coloc_per <- master.table.meta %>%
        as_tibble() %>%                     # turn the DFrame into a tibble
        group_by(module) %>%                # one group per module
        summarise(
            max_count = max(table(study_id)), # largest # of rows for any study in that module
            .groups = "drop"
        ) %>%
        filter(max_count != 1) %>%          # keep modules where that max ≠ 1
        pull(module)                        # extract just the module IDs

    max.module <- max(master.table.meta$module)

    for (module in modules_for_coloc_per) {

        abf.matrix <- flanders:::extract.ann.data.labf.matrix(
            sce.obj = ann.data,
            cs.list = master.table.meta$cs_name[master.table.meta$module == module]
        )
        res <- flanders:::coloc_percolation_matrix(
            finemap.data = abf.matrix,
            pp.H4.threshold = 0.5
        )

        if (!is.null(res$cp.object)) {
            for (cp_module in 1:length(res$cp.object$list.of.communities.labels)) {
                max.module <- max.module + 1
                cs.lab <- res$cp.object$list.of.communities.labels[[cp_module]]
                master.res <- colData(ann.data)[cs.lab, ]
                master.res$module <- max.module
                master.table.meta <- rbind(master.table.meta, master.res)
            }
        } else {
            max.module <- max.module + 1
            master.res <- colData(ann.data)[master.table.meta$cs_name[master.table.meta$module == module], ]
            master.res <- as.data.table(master.res)
            if (length(master.table.meta$cs_name[master.table.meta$module == module]) > 2) {
                master.res$pop <- apply(
                    t(master.res$cs_name), 2,
                    function(x) {
                        if (grepl("::", x)) {
                            unlist(strsplit(x, split = "::"))[1]
                        } else {
                            NA
                        }
                    }
                )
                master.res2 <- master.res %>%
                    group_by(pop) %>%
                    filter(
                        n() == 1 |
                        (pchisq((beta / se)^2, df = 1, lower.tail = FALSE, log.p = TRUE) ==
                            min(pchisq((beta / se)^2, df = 1, lower.tail = FALSE, log.p = TRUE)))
                    ) %>%
                    ungroup()
            } else {
                master.res2 = master.res
            }
            master.res2$module = max.module
            master.res2 = master.res2[, names(master.table.meta)]
            master.table.meta = rbind(master.table.meta, master.res2)
        }

    }

    # Remove rows from master.table.meta where the module is in partly_connected_modules using dplyr
    master.table.meta <- master.table.meta %>%
        filter(!(module %in% modules$partly_connected_modules$module))

    # Remove modules that have only one member using dplyr
    master.table.meta <- master.table.meta %>%
        group_by(module) %>%
        filter(n() > 1) %>%
        ungroup() %>% 
        as.data.table

    print(paste0("There are ", nrow(coloc_res_to_MA), " pairs of colocalized credible sets of same traits across studies."))

    # Convert each row to a list of two elements: hit1 and hit2
    cs_lists <- master.table.meta %>%
        group_by(module) %>%
        summarise(cs_names = list(cs_name), .groups = "drop") %>%
        pull(cs_names)

    i <- 0

    meta_results <- lapply(
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

    return(meta_results)
    
}