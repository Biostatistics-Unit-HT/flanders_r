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
4. [Additional Resources](#additional-resources)  
5. [Acknowledgments](#acknowledgments)

---

## Installation

### Simple Installation

To install the required R packages from CRAN, BioConductor and github, you can run the following commands in your R session:

```r
install.packages("data.table")
install.packages("dplyr")
install.packages("Matrix")
install.packages("optparse")

# You need these Bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("zellkonverter")
BiocManager::install("scrnaseq")

# Install the flanders package
# You need the devtools package if not installed.
install.packages("devtools")
devtools::install_github("anndata/anndataR")
devtools::install_github("Biostatistics-Unit-HT/flanders")
```
This guide installs all dependencies from CRAN and Bioconductor for a straightforward R-based setup.
### Installation via Conda Environment

For a more reproducible environment or if you need to interface with Python’s AnnData via reticulate, you can set up a Conda environment as follows:
```bash
mamba create -p <your_folder_with_conda_envs>/flanders_r \
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
  r-optparse=1.7.5

conda activate <your_folder_with_conda_envs>/flanders_r
git clone git@github.com:Biostatistics-Unit-HT/flanders_r.git
install.packages("flanders_r", repos = NULL, type = "source")

```

---

## Quick Start

### Scenario 1: Starting with an existing AnnData

If you already have an AnnData:
```r
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
---

## Function Reference

- **finemap2anndata**  
  Converts a set of finemapping `.rds` files into a single AnnData object with credible set metadata.

- **anndata2coloc_input**  
  Generates a data frame specifying pairs of credible sets for colocalization tests.

- **anndata2coloc**  
  Performs colocalization tests using the provided AnnData object and trait-pair table.

---

## Additional Resources

If the runtime for your colocalization tests becomes large, use the [flanders_nf_coloc](https://github.com/Biostatistics-Unit-HT/flanders_nf_coloc) Nextflow pipeline for scalable colocalization analysis.

---

## Acknowledgments

Contributions, bug reports, and feature requests are welcome. Open an issue in case of [issue](https://github.com/Biostatistics-Unit-HT/flanders_r/issues), bug or feature request.
