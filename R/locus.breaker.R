#' Identify significant association loci from GWAS summary statistics
#'
#' Groups variants passing a screening p-value threshold into loci based on
#' their chromosome and physical distance. Retains loci containing at least
#' one variant passing the significance threshold and reports the lead
#' variant for each retained locus.
#'
#' @param res A data.table containing GWAS summary statistics, with columns
#'   for chromosome, base-pair position, and p-value. Positions and p-values
#'   must be numeric. Additional columns are retained for the lead variant.
#' @param p.sig Numeric significance threshold. A locus is retained if its
#'   minimum p-value is strictly less than this value. Defaults to
#'   \code{5e-8}.
#' @param p.limit Numeric screening threshold. Only variants with p-values
#'   strictly less than this value are used to construct loci. Defaults to
#'   \code{1e-5}.
#' @param hole.size Numeric distance in base pairs. Consecutive variants
#'   passing \code{p.limit} on the same chromosome are assigned to separate
#'   loci when their distance is strictly greater than this value.
#'   Defaults to \code{250000}.
#' @param p.label Character string naming the p-value column.
#'   Defaults to \code{"P"}.
#' @param chr.label Character string naming the chromosome column.
#'   Defaults to \code{"CHR"}.
#' @param pos.label Character string naming the base-pair position column.
#'   Defaults to \code{"BP"}.
#'
#' @details
#' Variants are first sorted by chromosome and position, then filtered
#' using \code{p.limit}. Within each chromosome, gaps between consecutive
#' retained variants define locus boundaries. A gap equal to
#' \code{hole.size} does not split a locus.
#'
#' Each locus containing at least one variant with a p-value below
#' \code{p.sig} is reported. Its boundaries are the minimum and maximum
#' positions of the variants passing \code{p.limit} within that locus;
#' no flanking region is added. A locus containing a single variant has
#' identical start and end positions.
#'
#' The lead variant is the variant with the smallest p-value in the locus.
#' Ties are resolved by selecting the first variant in the sorted data.
#'
#' The input data.table is sorted by reference, which also changes the row
#' order of the object supplied by the caller. Use
#' \code{data.table::copy()} before calling this function if the original
#' row order must be preserved.
#'
#' @return A data.table with one row per significant locus, containing:
#'   \describe{
#'     \item{chr}{Chromosome of the locus.}
#'     \item{start}{Minimum position among the variants defining the locus.}
#'     \item{end}{Maximum position among the variants defining the locus.}
#'   }
#'   Remaining columns contain the lead variant's values from \code{res},
#'   excluding the original chromosome column specified by
#'   \code{chr.label}. The original position and p-value columns are retained.
#'
#'   Returns \code{NULL} if no significant loci are identified.
#'   Column types may be coerced during result assembly; for example,
#'   character annotations can cause numeric values to be returned as
#'   character strings.
#'
#' @import data.table
#' @export
locus.breaker <- function(
    res,
    p.sig = 5e-08,
    p.limit = 1e-05,
    hole.size = 250000,
    p.label = "P",
    chr.label = "CHR",
    pos.label = "BP"){
  
  #res <- as.data.frame(res)
  setorderv(res, c(chr.label, pos.label))
  res = res[get(p.label) < p.limit]
  trait.res = c()
  
  for(j in unique(res[, get(chr.label)])) {
    res.chr = res[get(chr.label) == j]
    if (nrow(res.chr) > 1) {
      pos_values <- res.chr[, get(pos.label)] # Get the pos.label column values
      holes <- pos_values[-1] - pos_values[-length(pos_values)]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr), 
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]), 
            ]
          }
          if (min(res.loc[, get(p.label)]) < p.sig) {
            start.pos = min(res.loc[, get(pos.label)], na.rm = T)
            end.pos = max(res.loc[, get(pos.label)], na.rm = T)
            chr = j
            best.snp = res.loc[which.min(res.loc[, get(p.label)]), 
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (min(res.loc[, get(p.label)]) < p.sig) {
          start.pos = min(res.loc[, get(pos.label)], na.rm = T)
          end.pos = max(res.loc[, get(pos.label)], na.rm = T)
          chr = j
          best.snp = res.loc[which.min(res.loc[, get(p.label)]), 
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (min(res.loc[, get(p.label)]) < p.sig) {
        start.pos = min(res.loc[, get(pos.label)], na.rm = T)
        end.pos = max(res.loc[, get(pos.label)], na.rm = T)
        chr = j
        best.snp = res.loc[which.min(res.loc[,get(p.label)]), 
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.table(trait.res, stringsAsFactors = FALSE)
    trait.res[, (chr.label) := NULL] # Remove the column specified by chr.label
    setnames(trait.res, 1:3, c("chr", "start", "end")) # Rename the first 3 columns
    setattr(trait.res, "row.names", NULL) # Remove row names
  }
  return(trait.res)
}
