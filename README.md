# pleisol

```r

reticulate::use_virtualenv("~/rstudio/virtualenvs/r-reticulate", required = TRUE)

anndata <- reticulate::import("anndata")

library(pleisol)
library(data.table)
library(dplyr)

# Define the folder containing the finemap files
finemap_folder <- "/group/pirastu/prj_014_huvec_coloc/2024_04_huvec_hcoloc_analysis_V2/results/finemap/"

# List all chr22 finemap files and filter out GWAS files
finemap_chr22_files <- list.files(finemap_folder, pattern = "chr22.*\\.rds", full.names = TRUE)
finemap_chr22_files <- finemap_chr22_files[!grepl("GWAS", finemap_chr22_files)]

# List all GWAS finemap files
finemap_gwas_files <- list.files(finemap_folder, pattern = "GWAS.*\\.rds", full.names = TRUE)

# Convert GWAS finemap files to AnnData format
gwas_ad <- finemap2anndata(
  finemap_files = finemap_gwas_files,
  output_file = "/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData/HUVEC_GWAS_combined_credible_sets.h5ad",
  panel = "HRC"
)

# Convert chr22 molecular QTL finemap files to AnnData format
chr22_molQTL_ad <- finemap2anndata(
  finemap_files = finemap_chr22_files,
  output_file = "/group/pirastu/prj_013_horizontal_codec/2024_06_20_AnnData/HUVEC_chr22_combined_credible_sets.h5ad",
  panel = "HRC"
)

chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")[414:416]

# Perform colocalization analysis on the GWAS AnnData object
gwas_ad.coloc <- anndata2coloc(gwas_ad)

# View summary of colocalization results
print(gwas_ad.coloc$summary)

# View detailed results by SNP
print(gwas_ad.coloc$results)
