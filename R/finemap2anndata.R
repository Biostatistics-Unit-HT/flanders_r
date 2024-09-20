#' Convert finemap files to AnnData
#'
#' This function reads finemap `.rds` files, filters rows where `is_cs` is `TRUE`,
#' and creates a sparse matrix of lABF values for SNPs across credible sets.
#' The function then converts the sparse matrix to an AnnData object and
#' includes metadata about the start and end positions of each credible set.
#'
#' @param finemap_files A character vector of file paths to the finemap `.rds` files.
#' @param study_id A character vector specifying the study_id for each of the finemap_files
#' @param phenotype_id A character vector specifying the phenotype_id for each of the finemap_files
#' @param output_file A character string specifying the path to save the resulting AnnData object as an `.h5ad` file.
#' @param panel A character string specifying the SNP genotyping/imputation panel
#'
#' @return The AnnData object containing the sparse matrix of lABF values and metadata.
#' @export
#' @import data.table
#' @import dplyr
#' @import anndata
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(anndata)
#' library(data.table)
#' library(dplyr)
#' library(Matrix)
#'
#' # Define the folder containing the finemap files
#' finemap_folder <- "/group/pirastu/prj_014_huvec_coloc/2024_04_huvec_hcoloc_analysis_V2/results/finemap/"
#'
#' # List all files matching the pattern
#' finemap_chr22_files <- list.files(finemap_folder, pattern = "chr22.*\\.rds", full.names = TRUE)
#' finemap_chr22_files <- finemap_chr22_files[!grepl("GWAS", finemap_chr22_files)]
#' # Function to extract study_id and phenotype_id from filename
#' get_study_phenotype_id <- function(file_name) {
#'   pattern <- "^(.*?chr[0-9XY]+)_(.*?)_chr[0-9XY]+:.*\\.rds$"
#'   matches <- regmatches(file_name, regexec(pattern, file_name))
#'   study_phenotype_id <- unlist(lapply(matches, function(x) x[2:3]))
#'   return(study_phenotype_id)
#' }
#'
#' # Collect study_id and phenotype_id for each credible set
#' study_phenotype_ids <- t(sapply(basename(finemap_chr22_files), get_study_phenotype_id))
#' colnames(study_phenotype_ids) <- c("study_id", "phenotype_id")
#' ad_chr_22 <- finemap2anndata(
#'   finemap_files = finemap_chr22_files,
#'   study_id = study_phenotype_ids$study_id,
#'   phenotype_id = study_phenotype_ids$phenotype_id,
#'   output_file = output_file_chr22,
#'   panel = "HRC"
#' )
#'
#' finemap_chr21_files <- list.files(finemap_folder, pattern = "chr21.*\\.rds", full.names = TRUE)
#' finemap_chr21_files <- finemap_chr21_files[!grepl("GWAS", finemap_chr21_files)]
#'
#' # List all ATAC files matching the pattern and filter out those containing "GWAS"
#' finemap_chr22_ATAC_files <- list.files(finemap_folder, pattern = "^ATAC.*\\chr22.*\\.rds", full.names = TRUE)
#' finemap_chr22_ATAC_files <- finemap_chr22_ATAC_files[!grepl("GWAS", finemap_chr22_ATAC_files)]
#'
#' # List all ChipSeq files matching the pattern and filter out those containing "GWAS"
#' finemap_chr22_ChipSeq_files <- list.files(finemap_folder, pattern = "^ChipSeq.*\\chr22.*\\.rds", full.names = TRUE)
#' finemap_chr22_ChipSeq_files <- finemap_chr22_ChipSeq_files[!grepl("GWAS", finemap_chr22_ChipSeq_files)]
#'
#' # Define the output file paths
#' output_file_chr22 <- file.path("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData", "HUVEC_chr22_combined_credible_sets.h5ad")
#' output_file_chr21 <- file.path("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData", "HUVEC_chr21_combined_credible_sets.h5ad")
#' output_file_chr22_ATAC <- file.path("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData", "HUVEC_chr22_ATAC_combined_credible_sets.h5ad")
#' output_file_chr22_ChipSeq <- file.path("/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData", "HUVEC_chr22_ChipSeq_combined_credible_sets.h5ad")
#'
#' ad_chr_22 <- finemap2anndata(
#'   finemap_files = finemap_chr22_files,
#'   output_file = output_file_chr22,
#'   panel = "HRC"
#' )
#'
#' ad_chr_21 <- finemap2anndata(
#'   finemap_files = finemap_chr21_files,
#'   output_file = output_file_chr21,
#'   panel = "HRC"
#' )
#'
#' ad_chr_22_ATAC <- finemap2anndata(
#'   finemap_files = finemap_chr22_ATAC_files,
#'   output_file = output_file_chr22_ATAC,
#'   panel = "HRC"
#' )
#'
#' ad_chr_22_ChipSeq <- finemap2anndata(
#'   finemap_files = finemap_chr22_ChipSeq_files,
#'   output_file = output_file_chr22_ChipSeq,
#'   panel = "HRC"
#' )
#'
#' chr_21_chr_22_outer <- concat(list(ad_chr_21, ad_chr_22), join = "outer")
#' chr_22_ATAC_chr_22_ChipSeq_outer <- concat(list(ad_chr_22_ATAC, ad_chr_22_ChipSeq), join = "outer")
#' }
finemap2anndata <- function(
    finemap_files,
    study_id = NULL,
    phenotype_id = NULL,
    output_file,
    panel = NULL
){

  # Initialize a list to store the filtered data
  filtered_data_list <- list()

  # Loop through each file and filter data
  print("Reading finemap files...")

  min_res_labf <- c()
  
  # Initialize counters and lists
  failed_files <- c()  # to store failed files
  success_count <- 0      # to count successful reads
  failed_count <- 0       # to count failed reads

  for (finemap_file in finemap_files) {
    
    # Try to read the .rds file and catch any errors or warnings
    result <- tryCatch({
      finemap <- data.table::data.table(readRDS(finemap_file))
      if("pC" %in% colnames(finemap))
        finemap <- finemap %>% rename(p=pC)
      
      success_count <- success_count + 1
      
      # If successful, continue with the rest of the operations
      min_res_labf <- c(min_res_labf, min(finemap$lABF))
      
      # Filter rows where is_cs is TRUE
      filtered_finemap <- finemap[is_cs == TRUE]
      
      # Append the filtered data to the list
      filtered_data_list[[basename(finemap_file)]] <- filtered_finemap
      
      NULL  # return NULL on success so tryCatch does not return a value
    }, error = function(e) {
      # On error, append the file name to the failed files list
      failed_files <<- c(failed_files,basename(finemap_file))
      failed_count <<- failed_count + 1
      #message(paste("Error in reading file:", basename(finemap_file), "- Skipping"))
      NULL  # return NULL on error
    })
  }

  
  study_id <- study_id[-sapply(failed_files,function(x) grep(x,finemap_files))]
  phenotype_id <- phenotype_id[-sapply(failed_files,function(x) grep(x,finemap_files))]
  
  # Output how many files were successfully read and how many failed
  cat("Number of successfully read files:", success_count, "\n")
  cat("Number of failed files:", failed_count, "\n")
  if (length(failed_files) > 0) {
    cat("Failed files:", paste(unlist(failed_files), collapse = ", "), "\n")
  } else {
    cat("Failed files: None\n")
  }

  # Collect all unique SNPs
  all_snps <- unique(unlist(lapply(filtered_data_list, function(x) x$snp)))

  element_indices <- stats::setNames(seq_along(all_snps), all_snps)

  # Collect all credible set names
  credible_sets <- names(filtered_data_list)

  # Create a sparse matrix to store lABF values
  lABF_matrix_sparse <- Matrix::Matrix(
    0,
    nrow = length(credible_sets),
    ncol = length(all_snps),
    sparse = TRUE,
    dimnames = list(credible_sets, all_snps)
  )

  print("Populating sparse matrix...")

  top_pvalue <- c()

  # Fill the sparse matrix with lABF values
  for (credible_set in credible_sets) {
    credible_data <- filtered_data_list[[credible_set]]
    row_index <- match(credible_set, credible_sets)
    col_indices <- element_indices[credible_data$snp]
    lABF_values <- credible_data$lABF
    lABF_matrix_sparse[row_index, col_indices] <- lABF_values
    top_pvalue <- c(top_pvalue,min(credible_data$p))
  }

  print("Creating AnnData object...")

  # Convert the sparse matrix to an AnnData object
  ad <- anndata::AnnData(X = lABF_matrix_sparse)

  print("Creating obs meta data...")

  # Extract start and end positions from the credible set names
  get_chr_start_end <- function(file_name) {
    pattern <- "locus_(chr[0-9XY]+)_([0-9]+)_([0-9]+)_susie_finemap.rds"
    matches <- regmatches(file_name, regexec(pattern, file_name))
    chr_start_end <- unlist(lapply(matches, function(x) x[2:4]))
    return(chr_start_end)
  }

  # Collect chromosome, start, and end positions for each credible set
  chr_start_end_positions <- t(sapply(credible_sets, get_chr_start_end))
  colnames(chr_start_end_positions) <- c("chr", "start", "end")

  # Convert to a data frame and set row names
  obs_df <- as.data.frame(chr_start_end_positions, stringsAsFactors = FALSE)
  obs_df$study_id <- study_id
  obs_df$phenotype_id <- phenotype_id
  obs_df$start <- as.numeric(obs_df$start)
  obs_df$end <- as.numeric(obs_df$end)
  obs_df$top_pvalue <- top_pvalue
  obs_df$min_res_labf <- min_res_labf
  obs_df$panel <- panel
  obs_df$cs_name <- credible_sets
  rownames(obs_df) <- credible_sets

  # Assign the data frame to ad$obs
  ad$obs <- obs_df

  #print("Writing AnnData to a disk...")

  # Save the AnnData object to .h5ad file
  #ad$write_h5ad(output_file)

  print("Done...")

  return(ad)
}
