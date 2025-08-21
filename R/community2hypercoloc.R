#' Perform hypercoloc analysis on a community of credible sets
#'
#' Extracts aligned beta and standard error matrices for a community of credible sets
#' and performs hypercoloc analysis to test for colocalization across multiple traits.
#'
#' @param sce_per_chr A list of SingleCellExperiment objects, one per chromosome.
#' @param community A list containing information about the community, including:
#'   \itemize{
#'     \item \code{nodes} - Character vector of credible set names
#'     \item \code{chr} - Chromosome number
#'   }
#'
#' @return A hypercoloc result object containing posterior probabilities for different
#' colocalization hypotheses across all traits in the community.
#'
#' @details The function first extracts the appropriate chromosome-specific SCE object,
#' then uses \code{extract_data_hypercoloc} to get aligned beta and SE matrices,
#' and finally performs hypercoloc analysis using the \code{hyprcoloc} package.
#'
#' @examples
#' \donttest{
#' library(SingleCellExperiment)
#' library(hyprcoloc)
#' library(data.table)
#' library(zellkonverter)
#' 
#' sce <- readH5AD(
#'   file = "/ssu/bsssu/anndata_finemapping_repository/cardinal_UKBB_europeans_ELGH_pakistani_cis_finemapping_full_FLIP_FIXED_20250611.h5ad",
#'   verbose = TRUE,
#'   reader = "R"
#' )
#' 
#' # Split SCE by chromosome
#' sce_by_chr <- pblapply(unique(colData(sce)$chr), function(chr) {
#'   snp_rows <- which(rowData(sce)$chr == chr)
#'   cred_cols <- which(colData(sce)$chr == chr)
#'   sce[snp_rows, cred_cols]
#' })
#' names(sce_by_chr) <- unique(colData(sce)$chr)
#' 
#' # Example community
#' community <- list(
#'   nodes = c("study1::trait1::cs1", "study2::trait1::cs2", "study3::trait1::cs3"),
#'   chr = 17
#' )
#' 
#' # Perform hypercoloc analysis
#' hyper_res <- community2hypercoloc(sce_by_chr, community)
#' }
#' }
#' @export
community2hypercoloc <- function(sce_per_chr, community) {
  # Subsetting sce with SNPs in list of hits
  sce_sm <- sce_per_chr[[as.character(community$chr)]]
  
  # Extract aligned beta and SE matrices
  dt_for_hprcoloc <- extract_data_hypercoloc(sce_sm, community)
  
  # Perform hypercoloc analysis
  hypr_res <- hyprcoloc::hyprcoloc(
    effect.est = dt_for_hprcoloc$beta_mat %>% as.matrix,
    effect.se = dt_for_hprcoloc$se_mat %>% as.matrix,
    trait.names = colnames(dt_for_hprcoloc$beta_mat),
    snp.id = row.names(dt_for_hprcoloc$beta_mat)
  )
  
  # ---- Build post_cp_tmp_module (cs_name, module) ----
  res <- hypr_res$results
  if (!is.data.frame(res)) stop("hypr_res$results not found or not a data.frame")
  
  rows_with_traits <- !is.na(res$traits) & res$traits != "None"
  rows_none        <- !is.na(res$traits) & res$traits == "None"
  
  pieces <- list()
  
  # Rows where traits != "None": split 'traits' by ", " and attach iteration as module
  if (any(rows_with_traits)) {
    for (k in which(rows_with_traits)) {
      cs_vec <- unlist(strsplit(res$traits[k], ",\\s*"))
      if (length(cs_vec)) {
        pieces[[length(pieces) + 1]] <- data.frame(
          cs_name = cs_vec,
          module  = res$iteration[k],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Rows where traits == "None": use dropped_trait (can be comma-separated too)
  if (any(rows_none)) {
    for (k in which(rows_none)) {
      drp <- res$dropped_trait[k]
      if (!is.na(drp) && nzchar(drp)) {
        cs_vec <- unlist(strsplit(drp, ",\\s*"))
        pieces[[length(pieces) + 1]] <- data.frame(
          cs_name = cs_vec,
          module  = res$iteration[k],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  post_cp_tmp_module <-
    if (length(pieces)){
      do.call(rbind, pieces)
      } else {
      data.frame(cs_name = character(0), module = integer(0), stringsAsFactors = FALSE)
    }
  
  # Replace module values with random strings while preserving groupings
  if (nrow(post_cp_tmp_module) > 0) {
    unique_modules <- unique(post_cp_tmp_module$module)
    module_mapping <- setNames(
      stringi::stri_rand_strings(
        n = length(unique_modules), 
        length = 30, 
        pattern = "[A-Za-z0-9]"
      ),
      unique_modules
    )
    post_cp_tmp_module$module <- module_mapping[as.character(post_cp_tmp_module$module)]
  }
  
  # ---- Build coloc.matrix (0/1, 1 on diagonal, 1 if in same module) ----
  if (nrow(post_cp_tmp_module) > 0) {
    cs_levels <- unique(post_cp_tmp_module$cs_name)
    N <- length(cs_levels)
    coloc_matrix <- matrix(0L, nrow = N, ncol = N,
                           dimnames = list(cs_levels, cs_levels))
    diag(coloc_matrix) <- 1L
    
    # set blocks to 1 for members in the same module
    idx_by_module <- split(seq_len(nrow(post_cp_tmp_module)), post_cp_tmp_module$module)
    for (idx in idx_by_module) {
      members <- post_cp_tmp_module$cs_name[idx]
      pos <- match(members, cs_levels)
      coloc_matrix[pos, pos] <- 1L
    }
  } else {
    coloc_matrix <- matrix(integer(0), nrow = 0, ncol = 0,
                           dimnames = list(character(0), character(0)))
  }
  
  # Create abf.list output
  labf_matrix <- extract.ann.data.labf.matrix(
    sce_sm, 
    community$nodes
  )
  
  labf_matrix <- apply(
    labf_matrix,
    2,
    icoloc.regularization
  )
  
  cl.list <- list(list(clique.map=NULL))
  cl.list[[1]]$clique.map <- split(
    f = post_cp_tmp_module$module,
    x = match(post_cp_tmp_module$cs_name, colnames(labf_matrix))
  )
  
  
  abf.list <- extract.abf.list.matrix(
    data = labf_matrix,
    clique.list = cl.list
  )
  
  # ---- Assemble required return shape ----
  inner <- list(
    abf.list           = abf.list,               # per spec: "let's forget about it for now"
    coloc.matrix       = coloc_matrix,
    post_cp_tmp_module = post_cp_tmp_module
  )
  
  # Function must return a list of single element whose only element
  # is this 3-field list:
  list(inner)
  
}