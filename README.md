# pleisol

```r

## Create conda environment with mamba

mamba create -p <your_folder_with_conda_envs>/pleisol -c conda-forge -c bioconda -c R \
r-base=4.3 r-susier=0.12.35 r-coloc=5.2.3 r-data.table=1.15.4 r-dplyr=1.1.4 r-anndata=0.7.5.6 r-mixsqp=0.3_54


For MacOS arm:
mamba create -p your_folder_with_conda_envs>/pleisol -c dnachun -c conda-forge -c bioconda -c R \
r-base=4.3 r-susier=0.12.35 r-coloc=5.2.3 r-data.table=1.15.4 r-dplyr=1.1.4 r-anndata=0.7.5.6 r-mixsqp=0.3_54

After creation of conda env, you should additionally install matrix package from 
https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.tar.gz

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

library(stringr)

chr22_molQTL_ad$obs$study_id <- str_extract(chr22_molQTL_ad$obs$cs_name, "^[A-Za-z]+_chr[0-9]+")

chr22_molQTL_ad$obs$phenotype_id <- str_match(chr22_molQTL_ad$obs$cs_name, "chr[0-9]+_([^_]+_[0-9]+|ENSG[0-9]+)")[,2]

# Perform colocalization analysis on the GWAS AnnData object
gwas_ad.coloc <- anndata2coloc(gwas_ad)

# View summary of colocalization results
print(gwas_ad.coloc$summary)

# View detailed results by SNP
print(gwas_ad.coloc$results)
```

## Data Column Specifications

### var (Variables)

| Column Name | Format / Content           | Description                                                                                                                                                         |
| ----------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `snp`       | `chr{num}:{pos}:EA:RA`    | **chr{num}**: chromosome identifier (any string)<br>**pos**: SNP position (in bp)<br>**EA**: effective allele (linked to beta)<br>**RA**: reference allele           |
| `chr`       | String (e.g. `"chr{num}"`)  | Chromosome where the SNP is located. Usually `"chr{num}"` format, but can be any string.                                                                            |
| `pos`       | Numeric (bp)               | Position of the SNP in base pairs (bp) on the physical map of the genome.                                                                                           |

**Note:** The row names of `ad$var` should be exactly equal to the values in `ad$var$snp`.

---

### obs (Observations)

#### Absolutely Necessary Columns

| Column Name    | Format / Content                            | Description                                                                                                                                                                                                                         |
| -------------- | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cs_name`      | `chr{num}::{study_id}::{trait_id}::{snp}`   | **chr{num}**: chromosome (formatted as in `var`)<br>**study_id**: identifier for the study<br>**trait_id**: refers to the trait (equivalent to `phenotype_id`)<br>**snp**: SNP with the highest logABF within the credible set (formatted as in `ad$var$snp`).            |
| `chr`          | String (e.g. `"chr{num}"`)                   | Chromosome where the credible set is located (same format as the `chr` field in `var`).                                                                                                                                             |
| `start`        | Numeric (bp)                                | Start position (in bp) of the analyzed locus, representing the beginning of the locus used for fine mapping.                                                                                                                        |
| `end`          | Numeric (bp)                                | End position (in bp) of the analyzed locus.                                                                                                                                                                                         |
| `study_id`     | String                                      | Identifier for the study.                                                                                                                                                                                                           |
| `phenotype_id` | String                                      | Identifier for the trait/phenotype analyzed within the corresponding study.                                                                                                                                                         |
| `min_res_labf` | Numeric (log-scale)                         | Minimal value of logABF in the locus. If logABF for all SNPs is not available, approximate using:<br><br>``[logsum(logABF)/coverage] - log(N_snps - N_CS_SNPs)``<br><br>- **coverage**: requested coverage (usually 99% or 95%)<br>- **logsum(logABF)**: log of the sum of ABF for SNPs in the credible set (using a log-sum to avoid overflow from extremely large values)<br>- Subtracting **log(N_snps - N_CS_SNPs)** gives the log(mean(ABF)) among SNPs outside the credible set. |

**Note:** The row names of `ad$obs` should be exactly equal to the values in `ad$obs$cs_name`.

---

#### Highly Advised to Have

| Column Name       | Format / Content     | Description                                                                                               |
| ----------------- | -------------------- | --------------------------------------------------------------------------------------------------------- |
| `min.abs.corr`    | Numeric              | Purity metric for the credible set: minimal absolute correlation between SNPs within the credible set.    |
| `mean.abs.corr`   | Numeric              | Mean absolute correlation among SNPs within the credible set.                                             |
| `median.abs.corr` | Numeric              | Median absolute correlation among SNPs within the credible set.                                           |

---

#### Good to Have

| Column Name   | Format / Content   | Description                                                                                            |
| ------------- | ------------------ | ------------------------------------------------------------------------------------------------------ |
| `top_pvalue`  | Numeric            | Lowest (either nominal or conditional) P-value of genetic association in the locus (useful for filtering). |
| `panel`       | String             | Name of the imputation panel used.                                                                     |
