% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anndata2coloc_row.R
\name{anndata2coloc_row}
\alias{anndata2coloc_row}
\title{Perform Pairwise Colocalization Analysis for a Specific Row}
\usage{
anndata2coloc_row(coloc_combo_row, ad_by_chr)
}
\arguments{
\item{coloc_combo_row}{A row from the colocalization guide table containing
information on trait pairs to be tested.}

\item{ad_by_chr}{A list of AnnData object, split by chromosome. Each element
of the list should have name like chr1, chr 21 and so on}
}
\value{
A list with colocalization results, including a summary with the top
SNPs and trait information.
}
\description{
This function performs colocalization analysis for a specific pair of traits
as defined in a row of the colocalization guide table. It retrieves relevant
SNPs within a specified region, imputes missing values, and calculates
colocalization probabilities between two datasets.
}
\examples{
\dontrun{
coloc_input <- anndata2coloc_input(ad)
coloc_combo_row <- coloc_input[1,]
# Example usage of anndata2coloc_row
result <- anndata2coloc_row(coloc_combo_row, ad_by_chr)
}
}
