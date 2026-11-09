#' Extract aligned beta and standard error matrices for hypercoloc analysis
#'
#' Extracts beta coefficients and standard errors from an AnnData or SingleCellExperiment object
#' for a specific community of credible sets, preparing matrices for hypercoloc analysis.
#'
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing summary statistics.
#' @param community A list containing information about the community, including:
#'   \itemize{
#'     \item \code{nodes} - Character vector of credible set names
#'     \item \code{chr} - Chromosome number
#'   }
#'
#' @return A list containing:
#'   \item{beta_mat}{Matrix of beta coefficients (SNPs x credible sets)}
#'   \item{se_mat}{Matrix of standard errors (SNPs x credible sets)}
#'
#' @details The function handles zero values by replacing them with small non-zero values
#' (0.00001 for beta, 1 for standard errors) to avoid numerical issues in downstream analysis.
#' It extracts SNPs in a region spanning from 100kb before the earliest start position to 
#' 100kb after the latest end position of all credible sets in the community.
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' library(data.table)
#' library(zellkonverter)
#' 
#' sce <- readH5AD(
#'   file = "/ssu/bsssu/anndata_finemapping_repository/cardinal_UKBB_europeans_ELGH_pakistani_cis_finemapping_full_FLIP_FIXED_20250611.h5ad",
#'   verbose = TRUE,
#'   reader = "R"
#' )
#' 
#' coloc_res <- fread("/ssu/bsssu/cardinal_finemap_coloc/anndata/anndata_020525_coloc/results/cardinal_UKBB_europeans_ELGH_pakistani_colocalization_table_all_merged.txt")
#' 
#' coloc_input <- anndata2coloc_input(sce)
#' 
#' # Create an edge list from your data
#' edges <- as.matrix(coloc_input[, c("t1", "t2")])
#' 
#' # Use both t1 and t2 because your nodes can appear in either
#' chr_lookup <- setNames(coloc_input$chr, coloc_input$t1)
#' chr_lookup2 <- setNames(coloc_input$chr, coloc_input$t2)
#' chr_lookup <- c(chr_lookup, chr_lookup2)
#' 
#' # Create the graph (undirected, since connections are bidirectional)
#' g <- igraph::graph_from_edgelist(edges, directed = FALSE)
#' 
#' # Find connected components (communities)
#' comps <- igraph::components(g)
#' 
#' # Get a list of lists, each being one community
#' community_list <- split(names(comps$membership), comps$membership)
#' 
#' community_chr <- sapply(community_list, function(members) chr_lookup[members[1]])
#' 
#' # As a list of lists:
#' community_info <- pbmapply(
#'   function(nodes, chr) list(nodes=nodes, chr=chr),
#'   community_list, community_chr,
#'   SIMPLIFY = FALSE
#' )
#' 
#' sce_by_chr <- pblapply(unique(coloc_input$chr), function(chr) {
#'   snp_rows <- which(rowData(sce)$chr == chr)
#'   cred_cols <- which(colData(sce)$chr == chr)
#'   sce[snp_rows, cred_cols]
#' })
#' 
#' names(sce_by_chr) <- unique(coloc_input$chr)
#' 
#' # Extract data for hypercoloc
#' data_matrices <- extract_data_hypercoloc(
#'   sce_by_chr[[as.character(community_info[[1]]$chr)]], 
#'   community_info[[1]]
#' )
#' }
#' @export
extract_data_hypercoloc <- function(ad_or_sce, community) {
  
  is_sce <- inherits(ad_or_sce, "SingleCellExperiment")
  
  # For simplicity, assume all credible sets are from the same chromosome
  chr_name <- as.character(community$chr)
  
  ad_obs <- if (is_sce) colData(ad_or_sce) else ad_or_sce$obs
  ad_var <- if (is_sce) rowData(ad_or_sce) else ad_or_sce$var
  
  # Extract start and end for each credible set
  starts <- ad_obs[community$nodes, "start"]
  ends   <- ad_obs[community$nodes, "end"]
  
  loc_start <- min(starts) - 1e5
  loc_end   <- max(ends) + 1e5
  
  # Find all SNPs in the region
  snp_pos   <- ad_var$pos
  snp_names <- rownames(ad_var)
  snps_for_coloc <- snp_names[snp_pos >= loc_start & snp_pos <= loc_end]
  
  # Extract the beta and SE matrices (SNPs x credible sets)
  if (is_sce) {
    rows <- match(snps_for_coloc, rownames(rowData(ad_or_sce)))
    cols <- match(community$nodes, rownames(colData(ad_or_sce)))
    if (any(is.na(rows))) {
      warning("Some SNPs missing in rowData; omitting.")
      valid <- !is.na(rows)
      snps_for_coloc <- snps_for_coloc[valid]
      rows <- rows[valid]
    }
    if (any(is.na(cols))) stop("Some credible set names not found in colData.")
    beta_mat <- assay(ad_or_sce, "beta")[rows, cols, drop = FALSE]
    se_mat <- assay(ad_or_sce, "se")[rows, cols, drop = FALSE]
  } else {
    beta_mat <- t(ad_or_sce[community$nodes, snps_for_coloc]$layers$beta)
    se_mat <- t(ad_or_sce[community$nodes, snps_for_coloc]$layers$se)
  }
  
  # Set row and column names
  colnames(beta_mat) <- community$nodes
  rownames(beta_mat) <- snps_for_coloc
  colnames(se_mat) <- community$nodes
  rownames(se_mat) <- snps_for_coloc
  
  # Handle zero values
  beta_mat[beta_mat == 0] <- 0.00001
  se_mat[se_mat == 0] <- 1
  
  return(
    list(
      beta_mat = beta_mat,
      se_mat = se_mat
    )
  )
}