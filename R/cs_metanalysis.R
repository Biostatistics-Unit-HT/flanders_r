#' Meta-analysis of credible sets from single-cell data
#'
#' Performs a meta-analysis on a list of credible sets (CS) from an AnnData or SingleCellExperiment object,
#' combining effect sizes and standard errors across studies or datasets.
#'
#' @param ad_or_sce An AnnData or SingleCellExperiment object containing summary statistics.
#' @param cs.list A character vector of credible set names to include in the meta-analysis.
#' @param coverage Numeric, the cumulative posterior probability coverage threshold (default: 0.99).
#'
#' @return A list containing:
#'   \item{cs.name}{A string with the meta-analysis credible set name.}
#'   \item{meta.results}{A data.frame with SNPID, meta-analyzed beta, se, logp, labf, and PPI.}
#'   \item{min.abf}{Minimum approximate Bayes factor in the credible set.}
#'   \item{start}{Start coordinate of the region.}
#'   \item{end}{End coordinate of the region.}
#'   \item{chr}{Chromosome of the region.}
#'   \item{coverage}{Sum of PPI for the credible set.}
#'   \item{logsum}{Sum of log Bayes factors for the credible set.}
#'   \item{best.snp.info}{A data.frame with the best SNP, its beta, and se.}
#'
#' @examples
#' \donotrun {
#' library(SingleCellExperiment)
#' library(data.table)
#' library(zellkonverter)
#' library(dplyr)
#'
#' ann.data <- readH5AD(
#'   file = "/ssu/bsssu/anndata_finemapping_repository/cardinal_UKBB_europeans_ELGH_pakistani_cis_finemapping_full.h5ad",
#'   verbose = TRUE,
#'   reader = "R"
#' )
#'
#' coloc_res <- fread("/ssu/bsssu/cardinal_finemap_coloc/anndata/anndata_020525_coloc/results/cardinal_UKBB_europeans_ELGH_pakistani_colocalization_table_all_merged.txt")
#'
#' # Select coloc results for same traits and different studies with high posterior probability of colocalization
#' coloc_res_to_MA <- coloc_res %>% filter(t1==t2 & t1_study_id != t2_study_id & PP.H4.abf >= 0.75)
#'
#' cs.metanalysis(
#'   ann.data,
#'   cs.list = c(
#'     "cardinal_ELGH_pakistani::likely_CD4_T_Naive:ENSG00000183172::chr22:42074313:C:T::L1",
#'     "cardinal_UKBB_europeans::likely_CD4_T_Naive:ENSG00000183172::chr22:42079699:G:T::L1"
#'   )
#' )
#'
#' metanalysis <- list()
#' for(i in 1:length(unique(master.table.meta$module))){
#'   mod.num <- unique(master.table.meta$module)[i]
#'   metanalysis[[i]] <- cs.metanalysis(
#'     ann.data,
#'     cs.list = master.table.meta$cs_name[master.table.meta$module == mod.num]
#'   )
#'   names(metanalysis)[i] <- paste0("Module_", mod.num)
#' }
#'}
#' @export

cs_metanalysis <- function(
  ad_or_sce = ann.data,
  cs.list = master.table.meta$cs_name[master.table.meta$module==1],
  coverage=0.99){

    print(cs.list)

    ad_sm <- slice_anndata(
        ad_or_sce = ad_or_sce,
        list_of_cs = cs.list
        )

    abf <- rowSums(ad_sm@assays@data$X)
    names(abf) <- rownames(ad_sm)

     PPI <- logbf_to_pp(abf)
     cs <- sort(PPI,decreasing=T)
     cs <- cs[1:(which(cumsum(cs)>=coverage)[1])]
     labfs <- abf[names(cs)]
     betas <- ad_sm[names(cs)]@assays@data$beta
     ses <- ad_sm[names(cs)]@assays@data$se
     weights <- 1/(ses^2)
     meta.beta <- (betas*weights)

    boundries <- c(
        start = min(ad_sm$start),
        end = max(ad_sm$end)
        )

    chr <- unique(ad_sm$chr)

     if(!is.null(nrow(meta.beta))){
       meta.beta <- rowSums(meta.beta)/rowSums(weights)
       meta.se <- sqrt(1/rowSums(weights))
       meta.p <- pchisq((meta.beta/meta.se)^2,df=1,lower=F,log.p =T )*-1
       meta.results <- data.frame(
        SNPID = names(cs),
        beta = meta.beta,
        se = meta.se,
        logp = meta.p,
        labf = abf[names(cs)],
        PPI = PPI[names(cs)]
        )
   
     }else{
       meta.beta <- sum(meta.beta) / sum(weights)
       meta.se <- sqrt(1/sum(weights))
       meta.p <- pchisq((meta.beta/meta.se)^2,df=1,lower=F,log.p =T )*-1
       meta.results <- data.frame(
        SNPID = names(cs),
        beta = meta.beta,
        se = meta.se,
        logp = meta.p,
        labf = abf[names(cs)],
        PPI = PPI[names(cs)]
        )
       
     }
     
     #cs.name=cs.data.list$locus.info$chr
     study.id <- paste(c("Meta",apply(t(cs.list),2,function(x)unlist(strsplit(x,split="::"))[1])),collapse=":")
     trait <- unlist(strsplit(cs.list[1],split="::"))[2]
     
     best.snp <- names(cs)[1]
     cs.name <- paste(c(study.id,trait,best.snp),collapse="::")
     best.snp.info <- cbind(
      best.snp,
      as.numeric(ad_sm[best.snp,]@assays@data$beta),
      as.numeric(ad_sm[best.snp,]@assays@data$se)
      )
       
       
     
     colnames(best.snp.info)=c("SNPID","beta","se")
     results = list(
        list_of_cs_metaanalyzed = cs.list,
        cs.name = cs.name,
        meta.results = meta.results,
        min.abf = min(abf),
        start = boundries["start"],
        end = boundries["end"],
        chr = chr,
        coverage = sum(PPI[names(labfs)]),
        logsum = sum(abf[names(labfs)]),
        best.snp.info = best.snp.info
    )

    results

}