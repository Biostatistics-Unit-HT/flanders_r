# AI Coding Instructions for Flanders

## Project Overview
**Flanders** is a statistical genetics R package for fine-mapping, colocalization, and pleiotropy analysis. It converts fine-mapping outputs into AnnData objects and performs pairwise colocalization tests to identify modules of co-regulated genetic signals.

## Key Architecture Concepts

### Core Data Flow
1. **Input**: `.rds` files from nf-flanders finemapping pipeline containing log Bayes factors (lABF)
2. **Conversion**: `finemap2anndata()` creates sparse matrices in AnnData format  
3. **Analysis**: Three-step colocalization workflow ending in joint fine-mapping

### AnnData Structure (Critical for AI Understanding)
- **Rows (var)**: SNPs with columns `snp`, `chr`, `pos` following format `chr{num}:{pos}:EA:RA`
- **Columns (obs)**: Credible sets with `cs_name` format `chr{num}::{study_id}::{trait_id}::{snp}`
- **Main assay X**: Sparse matrix of log approximated Bayes factors (lABF)
- **Additional layers**: `beta` (effect sizes), `se` (standard errors)

### Three-Phase Analysis Workflow
1. **Pairwise Colocalization**: `anndata2coloc_input()` → `anndata2coloc()` 
2. **Module Detection**: `identify_simple_modules()` uses igraph for community detection
3. **Joint Fine-mapping**: `run_joint_finemapping()` aggregates lABFs across modules

## Development Patterns

### R Package Structure
- Standard R package layout with `man/`, `R/`, `DESCRIPTION`, `NAMESPACE`
- Uses roxygen2 for documentation - run `make roxygen` before building
- Dependencies: Matrix (sparse), data.table, igraph, SingleCellExperiment, zellkonverter

### Critical Naming Conventions
- **Credible set names**: Follow `chr{num}::{study_id}::{trait_id}::{snp}` format exactly
- **SNP identifiers**: Use `chr{num}:{pos}:EA:RA` format consistently  
- **File patterns**: Fine-mapping files follow `locus_{chr}_{start}_{end}_finemap.rds`

### Sparse Matrix Operations
- Heavy use of `Matrix::dgCMatrix` for memory efficiency
- Pattern: `X@x`, `X@p`, `X@i` for direct sparse matrix element access
- Always use `Matrix::t()` for transposition of sparse matrices

### Error Handling Patterns
```r
# Standard pattern in finemap2anndata()
result <- tryCatch({
  finemap <- readRDS(finemap_file)
  # processing...
  NULL  # return NULL on success
}, error = function(e) {
  failed_files <<- c(failed_files, basename(finemap_file))
  NULL
})
```

## Key Functions to Understand

### Data Conversion
- `finemap2anndata()`: Main conversion function, handles file parsing and sparse matrix creation
- `fix_ad_var()`: Standardizes variant metadata across AnnData/SCE objects  
- `validate_anndata()`: Comprehensive validation with integration tests

### Colocalization Pipeline  
- `anndata2coloc_input()`: Matrix multiplication to find SNP overlaps between credible sets
- `anndata2coloc_row()`: Single pairwise colocalization using `icoloc.ht()`
- `identify_simple_modules()`: igraph-based community detection with connectivity checks

### Advanced Analysis
- `coloc_percolation_matrix()`: Refines partially-connected modules using recursive percolation
- `joint.finemapping()`: Aggregates lABFs across credible sets in a module
- `run_joint_finemapping()`: Batch wrapper for joint fine-mapping across all modules

## Build and Development

### Building Package
```bash
make roxygen  # Generate documentation
make build    # Create .tar.gz
make check    # R CMD check
make install  # Install locally
```

### Testing Integration
- Use `validate_anndata()` for comprehensive smoke tests
- Check against both SingleCellExperiment and Python AnnData formats
- Verify sparse matrix operations preserve memory efficiency

## Common Pitfalls to Avoid

1. **Matrix Orientation**: AnnData convention is samples×features, but genetic data is often transposed
2. **Sparse Matrix Operations**: Always use Matrix package functions, never convert to dense
3. **CS Name Parsing**: Regular expressions in `get_chr_start_end()` are format-specific
4. **Memory Management**: Process chromosome-wise for large datasets to avoid memory issues
5. **Module Numbering**: Ensure numeric module IDs don't collide when merging results

## External Dependencies
- **nf-flanders pipeline**: Upstream finemapping pipeline producing input `.rds` files
- **coloc package**: Core colocalization methods (Bayesian approach)
- **Python AnnData**: Cross-language compatibility via reticulate/zellkonverter
- **HyperCoLoc**: For multi-trait colocalization (github.com/jrs95/hyprcoloc)

When modifying this codebase, maintain the sparse matrix efficiency, preserve AnnData schema compatibility, and ensure chromosome-wise processing patterns for scalability.