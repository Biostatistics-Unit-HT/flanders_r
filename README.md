# flanders

**flanders** is an R package designed to seamlessly convert finemapping output files (e.g., from the [nf-flanders](https://github.com/Biostatistics-Unit-HT/Flanders) pipeline) into a unified AnnData object and facilitate colocalization analysis. The package provides functions to:

- Convert multiple  `*finemap.rds` files into a single AnnData object with credible set metadata.
- Generate an input table (`coloc_input`) for colocalization testing.
- Run pairwise colocalization tests, with minimal runtime overhead (typically 5–10 tests per second on standard hardware).

When processing small to moderate datasets, you can run colocalization tests on your PC or laptop. For large-scale analyses, consider using the [flanders_nf_coloc](https://github.com/Biostatistics-Unit-HT/flanders_nf_coloc) Nextflow pipeline.

---

## Table of Contents

1. [Installation](#installation)  
2. [Quick Start](#quick-start)  
   1. [Scenario 1: Starting with nf-flanders Finemapping Output](#scenario-1-starting-with-nf-flanders-finemapping-output)  
   2. [Scenario 2: Starting with an existing AnnData Object](#scenario-2-starting-with-an-existing-anndata-object)  
3. [Function Reference](#function-reference)  
4. [Additional Resources](#additional-resources)  
5. [Acknowledgments](#acknowledgments)

---

## Installation

Below are the steps to set up a conda environment with the required dependencies:


**For Linux/Windows:**
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
```

---

## Quick Start

### Scenario 1: Starting with nf-flanders Finemapping Output

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
    ad$write_h5ad("/path/to/output/credible_sets.h5ad")
```

3. **Generate Coloc Input Table. Write it if you want further run the nf-hcoloc pipeline**
```r
    coloc_input <- anndata2coloc_input(ad)
    fwrite(coloc_input, file = "/path/to/coloc_guide_table.csv")
```

4. **Perform Colocalization Analysis**
   ```r

    coloc_results <- anndata2coloc(ad, coloc_input)
    print(coloc_results)
   ```

### Scenario 2: Starting with an existing AnnData Object

If you already have an AnnData object:
```r
    coloc_input <- anndata2coloc_input(existing_anndata)
    coloc_results <- anndata2coloc(existing_anndata, coloc_input)
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
