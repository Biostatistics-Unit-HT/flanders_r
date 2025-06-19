#' Joint Fine-Mapping Across a Module of Credible Sets
#'
#' @description
#' This function performs joint fine-mapping across a group of credible sets (CSs) 
#' defined in a module. It aggregates log-approximated Bayes factors (lABFs) across 
#' all CSs, computes posterior inclusion probabilities (PIPs), and determines the 
#' 99% credible set based on cumulative PIP.
#'
#' @param ann.data.obj A `SingleCellExperiment` or similar annotated data object
#'        containing SNP-level lABFs (`X` assay), effect sizes (`beta` assay),
#'        and standard errors (`se` assay).
#' @param module.master.table A data.frame or data.table containing credible sets
#'        (CSs) for a given module. Must include a `cs_name` column that matches
#'        column names in `ann.data.obj`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{most.likely.snp.info}{A data.frame containing the most likely causal SNP
#'         (i.e., the top PIP SNP) with effect size (`beta`) and standard error (`se`)
#'         across all CSs in the module.}
#'   \item{cs.info}{A data.frame with SNPs in the 99% credible set including their
#'         posterior inclusion probabilities (`PIP`) and joint lABF scores.}
#'   \item{mim.labf}{The minimum lABF among SNPs not included in the 99% credible set.}
#'   \item{sum.labf}{The log-sum of all lABFs (used for normalizing posterior probabilities).}
#' }
#'
#' @details
#' This function first extracts the lABF matrix across all CSs in a given module using
#' `extract.ann.data.labf.matrix()`. It then computes joint evidence for each SNP by summing
#' lABFs across CSs. The SNPs are ranked by PIP and a 99% credible set is identified.
#'
#' The most likely causal SNP (i.e., top-ranked by PIP) is also returned with its associated
#' effect sizes and standard errors across CSs in the module.
#'
#' @seealso 
#' \code{\link{extract.ann.data.labf.matrix}}, \code{\link{logbf_to_pp}}, \code{\link[coloc.peRculation]{logsum}}
#'
#' @examples
#' \dontrun{
#' results <- joint.finemapping(
#'   ann.data.obj = ann.data,
#'   module.master.table = master_table[module == 3]
#' )
#' results$cs.info
#' }
#'
#' @importFrom data.table data.table
#' @export
#' 
joint.finemapping = function(ann.data.obj=ann.data,
                           module.master.table=NULL){
  
  
  abf.matrix <- extract.ann.data.labf.matrix(
    sce.obj = ann.data.obj,
    cs.list = module.master.table$cs_name)
  
  joint.fm <- rowSums(abf.matrix)
  PIP <- logbf_to_pp(joint.fm)
  
  cs.99 <- cumsum(sort(PIP,decreasing = T))
  cs.99 <- names(cs.99[1:(which(cs.99>0.99)[1])])
  
  most.likely.snp <- cs.99[1]
  betas <- ann.data.obj[most.likely.snp,module.master.table$cs_name]@assays@data$beta
  se <- ann.data.obj[most.likely.snp,module.master.table$cs_name]@assays@data$se
  
  most.likely.snp.info <- data.frame(
    cs_name = module.master.table$cs_name,
    most_likely_snp = most.likely.snp,
    beta = as.vector(betas),
    se = as.vector(se)
  )
  
  cs.info <- data.frame(
    snp.name = cs.99,
    PIP = PIP[cs.99],
    lABF = joint.fm[cs.99]
  )
  
  mim.labf <- min(joint.fm[which(!(joint.fm%in%cs.99))],na.rm=T)
  sum.labf <- coloc.peRculation:::logsum(joint.fm) 
  
  return(list(
    most.likely.snp.info = most.likely.snp.info,
    cs.info = cs.info,
    mim.labf = mim.labf,
    sum.labf = sum.labf)
  )
  
}
