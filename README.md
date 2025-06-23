# flanders

**flanders** is an R package designed to seamlessly convert finemapping output files (e.g., from the [nf-flanders](https://github.com/Biostatistics-Unit-HT/Flanders) pipeline) into a unified AnnData object and facilitate colocalization analysis. The package provides functions to:

- Convert multiple  `*finemap.rds` files into a single AnnData object with credible set metadata.
- Generate an input table (`coloc_input`) for colocalization testing.
- Run pairwise colocalization tests, with minimal runtime overhead (typically 5–10 tests per second on standard hardware).

When processing small to moderate datasets, you can run colocalization tests on your PC or laptop. For large-scale analyses, consider using the [flanders_nf_coloc](https://github.com/Biostatistics-Unit-HT/flanders_nf_coloc) Nextflow pipeline.

---

## Table of Contents

1. [Installation](#installation)  
   1. [Simple Installation](#simple-installation)  
   2. [Installation via Conda Environment](#installation-via-conda-environment)
2. [Quick Start](#quick-start)  
   1. [Scenario 1: Starting with an existing AnnData](#scenario-1-starting-with-an-existing-anndata)  
   2. [Scenario 2: Starting with nf-flanders Finemapping Output](#scenario-2-starting-with-nf-flanders-finemapping-output)   
3. [Function Reference](#function-reference)
4. [AnnData Column Specifications](#anndata-column-specifications)  
5. [Additional Resources](#additional-resources)  
6. [Acknowledgments](#acknowledgments)

---

## Installation

### Simple Installation

To install the required R packages from CRAN, BioConductor and github, you can run the following commands in your R session:

```r
install.packages("data.table")
install.packages("dplyr")
install.packages("Matrix")
install.packages("optparse")
install.packages("R.utils")
install.packages("stringi")
install.packages("qgraph")
install.packages("igraph")

# You need these Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("zellkonverter")
BiocManager::install("scRNAseq")

# Install the flanders package
# You need the devtools package if not installed.
install.packages("devtools")
install.packages("anndata")
devtools::install_github("Biostatistics-Unit-HT/flanders_r")
```
This guide installs all dependencies from CRAN and Bioconductor for a straightforward R-based setup.
### Installation via Conda Environment

For a more reproducible environment or if you need to interface with Python’s AnnData via reticulate, you can set up a Conda environment as follows:
```bash
mamba create -p /ssu/bsssu/conda_envs_bsssu/flanders_r_2.0 \
  -c conda-forge -c bioconda -c R \
  r-base=4.4 \
  bioconductor-singlecellexperiment=1.28.0 \
  bioconductor-zellkonverter=1.16.0 \
  bioconductor-scrnaseq=2.20.0 \
  r-susier=0.12.35 \
  r-coloc=5.2.3 \
  r-data.table=1.17.0 \
  r-dplyr=1.1.4 \
  r-anndata=0.7.5.6 \
  r-matrix=1.7_3 \
  r-optparse=1.7.5 \
  r-roxygen2=7.3.2 \
  r-igraph=2.0.3 \
  r-qgraph=1.9.8 \
  r-r.utils=2.13.0 \
  r-stringi

conda activate <your_folder_with_conda_envs>/flanders_r
git clone git@github.com:Biostatistics-Unit-HT/flanders_r.git
install.packages("flanders_r", repos = NULL, type = "source")

```

---

## Quick Start

### Scenario 1: Starting with an existing AnnData

If you already have an AnnData:
```r
   library(zellkonverter)
   library(SingleCellExperiment)
   library(scRNAseq)

   sce <- readH5AD("/path/to/output/my_anndata.h5ad",reader="R")
   coloc_input <- anndata2coloc_input(sce)
   coloc_results <- anndata2coloc(sce, coloc_input)

   print(coloc_results)
```

If you already have multiple AnnDatas:
```r
   list_of_ads <- list("/path/to/ad1.h5ad","/path/to/output/ad2.h5ad")
   list_of_ads <- lapply(list_of_ads,read_h5ad)

   merged_ad <- concat(
     list_of_ads,
     join="outer",
     merge="first"
   )
   merged_ad <- fix_ad_var(merged_ad)
   write_h5ad(merged_ad,"/path/to/merged_ad.h5ad")

   library(zellkonverter)
   library(SingleCellExperiment)
   library(scRNAseq)

   sce <- readH5AD("/path/to/merged_ad.h5ad",reader="R")
   coloc_input <- anndata2coloc_input(sce)
   coloc_results <- anndata2coloc(sce, coloc_input)
   print(coloc_results)
```

### Scenario 2: Starting with nf-flanders Finemapping Output

If you do not have an AnnData object yet:

1. **Convert Finemapping Files to AnnData**
```r

    library(flanders)
    library(data.table)
    library(dplyr)
    
    finemap_folder <- "/path/to/finemap/results/"
    finemap_files <- list.files(finemap_folder, pattern = "*\\.rds", full.names = TRUE)

    # create a vector of phenotype ids. Each element corresponds to each file of finemap_files.
    # phenotype_id <-

    # create a vector of study ids. Each element corresponds to each file of finemap_files.
    # study_id <- 
    
    ad <- finemap2anndata(
      finemap_files = finemap_files,
      phenotype_id = phenotype_id,
      study_id = study_id
    )
    
    # Optionally write the AnnData object to disk
    ad$write_h5ad("/path/to/output/my_anndata.h5ad")
```

3. **Generate Coloc Input Table. Write it if you want further run the nf-hcoloc pipeline**
```r
   library(zellkonverter)
   library(SingleCellExperiment)
   library(scRNAseq)
   sce <- readH5AD("/path/to/output/my_anndata.h5ad",reader="R")
   coloc_input <- anndata2coloc_input(sce)
   fwrite(coloc_input, file = "/path/to/coloc_guide_table.csv")
```

4. **Perform Colocalization Analysis**
```r
   coloc_results <- anndata2coloc(sce, coloc_input)
   print(coloc_results)
```

5. **Identify simple modules/group of co-regulated signals and create credible sets-wise master table**
```r
   library(igraph)

   modules <- identify_simple_modules(coloc_results, pp.H4.threshold=0.75, check_connectivity=FALSE)
   master_table <- create_master_table(ad, simple_modules=modules)
   print(master_table)
   ```  

6. **Identify modules/group of co-regulated signals, use coloc percolation to dientagloe the most comple ones, create credible sets-wise master table and perform joint fine-mapping**
```r
   library(igraph)
   library(qgraph)
   library(stringi)

   # Identify fully connected and partially connected modules
   modules <- identify_simple_modules(coloc_results, pp.H4.threshold=0.75, check_connectivity=TRUE)
   
   # Run partially connected module through coloc percolation to refine them
   groups <- run_coloc_percolation(
     modules$partly_connected_modules,
     ad,
     pp.H4.threshold=0.75,
     icoloc.reg=TRUE)

   # Create master table
   master_table <- create_master_table(
     ad,
     post_cp_list = groups,
     fully_connected_modules = modules$fully_connected_modules
  )

   # Perform joint fine-mapping
   joint_finemapping_obj <- run_joint_finemapping(ad, master_table)

   joint_finemapping_obj$master_table # Access master table with modules, most likely SNP from joint fine-mapping
   joint_finemapping_obj$joint.lABF.matrix # Access sparse matrix collecting joint fine-mapping lABFs (per credible set and SNP)

   ```     
---

## Function Reference

- **finemap2anndata**  
  Converts a set of finemapping `.rds` files into a single AnnData object with credible set metadata.

- **anndata2coloc_input**  
  Generates a data frame specifying pairs of credible sets for colocalization tests.

- **anndata2coloc**  
  Performs colocalization tests using the provided AnnData object and trait-pair table.

- **identify_simple_modules**
  Defines modules (communities) of colocalizing signals using a graph-based approach. Optionally, checks whether each module is fully connected.

- **create_master_table**
  Constructs a master table where each credible set has been assigned to a module. It can be run directly on simple module or modules which have been further refined by coloc percolation.

- **run_coloc_percolation**
  Runs colocalization-guided percolation analysis to refine modules having credible sets which did not all colocalise among themselves in pair-wise tests.

- **run_joint_finemapping**
  Once modules have been defined, runs a joint fine-mapping to identify a single credible set for each module.

---

## AnnData Column Specifications

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
| `min_res_labf` | Numeric (log-scale)                         | Minimal value of logABF in the locus. If logABF for all SNPs is not available, approximate using:<br><br>``[logsum(logABF)/coverage] - log(N_snps - N_CS_SNPs)``
| coverage | Numeric | requested coverage (usually 99% or 95%)
| logsum(logABF)| Numeric| log of the sum of ABF for SNPs in the credible set (using a log-sum to avoid overflow from extremely large values) | Subtracting **log(N_snps - N_CS_SNPs)** gives the log(mean(ABF)) among SNPs outside the credible set. |

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

## Additional Resources

If the runtime for your colocalization tests becomes large, use the [flanders_nf_coloc](https://github.com/Biostatistics-Unit-HT/flanders_nf_coloc) Nextflow pipeline for scalable colocalization analysis.

---

## Acknowledgments

Contributions, bug reports, and feature requests are welcome. Open an issue in case of [issue](https://github.com/Biostatistics-Unit-HT/flanders_r/issues), bug or feature request.
