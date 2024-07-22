#' Perform Colocalization Analysis on AnnData Object
#'
#' This function performs colocalization analysis on an AnnData object. It computes the product of 
#' the sparse matrix with its transpose to identify shared elements, retrieves trait information, 
#' and performs colocalization testing on the identified pairs.
#'
#' @param ad An AnnData object containing genetic data.
#'
#' @return A list of colocalization results for each pair of traits.
#' @export
#'
#' @examples
#' \dontrun{
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#'
#' ad <- read_h5ad("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData/HUVEC_chr22_combined_credible_sets.h5ad")
#' ad.coloc <- anndata2coloc(ad[1:10,])
#'
#' # Store ALL the summary output in a data frame, adding tested traits column and SAVE 
#' only_summary_df <- as.data.frame(rbindlist(lapply(ad.coloc, function(x) { x$summary })))
#' }
anndata2coloc <- function(ad) {
  
  matrix_product <- ad$X %*% t(ad$X)
  
  # Set the diagonal to zero (as each vector will have 1s on diagonal)
  diag(matrix_product) <- 0
  
  shared_elements <- which(matrix_product > 0, arr.ind = TRUE)
  shared_elements_df <- as.data.frame(shared_elements)
  
  # Map row and column indices to credible set names
  shared_elements_df$row <- rownames(shared_elements_df)
  shared_elements_df$col <- colnames(matrix_product)[shared_elements_df$col]
  
  # Keep unique pairs
  shared_elements_unique <- unique(t(apply(shared_elements_df, 1, sort)))
  colnames(shared_elements_unique) <- c("t1", "t2")
  rownames(shared_elements_unique) <- NULL
  
  # Retrieve trait info
  coloc_combo <- merge(
    shared_elements_unique,
    ad$obs %>% select(-panel),
    by.x = "t1",
    by.y = "cs_name"
  )
  
  coloc_combo <- merge(
    coloc_combo,
    ad$obs %>% select(-panel),
    by.x = "t2",
    by.y = "cs_name",
    suffixes = c("_t1", "_t2")
  )
  
  # Remove pair testing different conditional dataset for the same trait (study_id + phenotype_id)
  coloc_combo <- coloc_combo %>%
    dplyr::filter(study_id_t1 != study_id_t2 | (study_id_t1 == study_id_t2 & phenotype_id_t1 != phenotype_id_t2)) %>%
    dplyr::select(-chr_t1, -chr_t2, -contains("start"), -contains("end"))
  
  colnames(coloc_combo) <- c(
    "t2", "t1", "t1_study_id", "t1_phenotype_id", "t1_top_pvalue",
    "t2_study_id", "t2_phenotype_id", "t2_top_pvalue"
  )
  
  ###### COLOCALIZATION ######
  
  # Split the dataframe into a list of rows
  coloc_combo_ls <- split(coloc_combo, seq(nrow(coloc_combo)))
  
  # Perform coloc!
  coloc.full <- lapply(coloc_combo_ls, function(coloc_combo_row) {
    
    # Load-in precomputed lABF  
    dataset1 <- ad$X[coloc_combo_row$t1,]
    dataset1 <- dataset1[which(dataset1 != 0)]
    
    dataset2 <- ad$X[coloc_combo_row$t2,]
    dataset2 <- dataset2[which(dataset2 != 0)]
    
    dataset1 <- impute.labf(
      snp = union(names(dataset1), names(dataset2)),
      cred.set = names(dataset1),
      cred.set.labf = dataset1,
      imputed.labf = ad$obs[coloc_combo_row$t1, "min_res_labf"]
    )
    
    dataset2 <- impute.labf(
      snp = union(names(dataset2), names(dataset2)),
      cred.set = names(dataset2),
      cred.set.labf = dataset2,
      imputed.labf = ad$obs[coloc_combo_row$t2, "min_res_labf"]
    )
    
    icoloc.res <- icoloc.ht(
      dataset1,
      dataset2
    )
    
    # Retrieve important info from file name
    t1 <- ifelse(is.na(coloc_combo_row$t1_phenotype_id), coloc_combo_row$t1_study_id, paste0(coloc_combo_row$t1_study_id, "_", coloc_combo_row$t1_phenotype_id))
    cojo_snp1 <- gsub(paste0(t1, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t1)
    
    t2 <- ifelse(is.na(coloc_combo_row$t2_phenotype_id), coloc_combo_row$t2_study_id, paste0(coloc_combo_row$t2_study_id, "_", coloc_combo_row$t2_phenotype_id))
    cojo_snp2 <- gsub(paste0(t2, "_(.*)_locus_.*_finemap.rds"), "\\1", coloc_combo_row$t2)
    
    # Add top cojo SNPs and traits
    icoloc.res$summary <- icoloc.res$summary %>%
      mutate(
        t1_study_id = coloc_combo_row$t1_sty_id, 
        t1 = t1, 
        t2_study_id = coloc_combo_row$t2_study_id,
        t2 = t2,
        hit1 = cojo_snp1,
        hit2 = cojo_snp2
      )
    
    icoloc.res
    
  })
  
  return(coloc.full)
}
