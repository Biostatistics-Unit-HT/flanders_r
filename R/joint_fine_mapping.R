#' Joint Fine-Mapping Across a Module of Credible Sets
#'
#' @description
#' Performs joint fine-mapping across a group of credible sets (CSs) within a module.
#' Aggregates log-approximated Bayes factors (lABFs) across CSs, computes posterior inclusion
#' probabilities (PIPs), and identifies a 99% credible set based on cumulative PIP.
#'
#' @param ann.data.obj A `SingleCellExperiment` or similar object containing SNP-level assays,
#'        including `X` (lABFs), `beta` (effect sizes), and `se` (standard errors).
#'
#' @param module.master.table A `data.frame` or `data.table` containing credible sets for
#'        a single module. Must include a `cs_name` column with names matching columns
#'        in `ann.data.obj`.
#'
#' @param full_labf_matrix (Optional) A precomputed matrix of lABFs across all CSs
#'        (columns = CSs, rows = SNPs). If provided, the function will subset this matrix
#'        for the current module instead of calling `extract.ann.data.labf.matrix()` internally.
#'        This is recommended for batch fine-mapping across multiple modules to improve performance.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item most.likely.snp.info: A data.frame with the top PIP SNP and its associated beta and se across all CSs in the module.
#'   \item cs.info: A data.frame of SNPs in the 99% credible set, with their PIPs and joint lABFs.
#'   \item mim.labf: The minimum lABF among SNPs not in the 99% credible set.
#'   \item sum.labf: The log-sum of all lABFs, used to normalize posterior probabilities.
#' }
#' 
#' @details
#' If `full_labf_matrix` is not provided, the function calls `extract.ann.data.labf.matrix()`
#' on the current module’s CSs. To avoid redundant extraction in batch runs (e.g. across many modules),
#' pass in `full_labf_matrix` and let the function subset the relevant CSs internally.
#'
#' The top-ranked SNP by PIP is used to extract effect sizes and standard errors from the `beta`
#' and `se` assays of `ann.data.obj`.
#'
#' @seealso 
#' \code{\link{extract.ann.data.labf.matrix}}, 
#' \code{\link{logbf_to_pp}}, 
#' \code{\link[coloc]{logsum}}
#'
#' @examples
#' \dontrun{
#' # Run fine-mapping on module 3
#' results <- joint.finemapping(
#'   ann.data.obj = ann.data,
#'   module.master.table = master_table[module == 3]
#' )
#'
#' # Provide precomputed matrix for efficiency
#' labf_mat <- extract.ann.data.labf.matrix(ann.data, master_table$cs_name)
#' results <- joint.finemapping(
#'   ann.data.obj = ann.data,
#'   module.master.table = master_table[module == 3],
#'   full_labf_matrix = labf_mat
#' )
#' }
#'
#' @importFrom data.table data.table
#' @export
#' 
joint.finemapping = function(ann.data.obj,
                             module.master.table,
                             full_labf_matrix = NULL){
  
  # Subset the full_labf_matrix to the cs_names for this module
  if (is.null(full_labf_matrix)) {
    abf.matrix <- extract.ann.data.labf.matrix(
      sce.obj = ann.data.obj,
      cs.list = module.master.table$cs_name
    )
  } else {
    abf.matrix <- full_labf_matrix[, module.master.table$cs_name, drop = FALSE]
  }
  
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
  sum.labf <- coloc:::logsum(joint.fm) 
  
  return(list(
    most.likely.snp.info = most.likely.snp.info,
    cs.info = cs.info,
    mim.labf = mim.labf,
    sum.labf = sum.labf)
  )
  
}
