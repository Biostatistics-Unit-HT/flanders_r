slice_anndata <- function(
    ad_or_sce = ad_or_sce,
    list_of_cs = list_of_cs)
    {

  ad_or_sce <- ad_or_sce[,list_of_cs]

  boundries <- c(
    start = min(ad_or_sce$start),
    end = max(ad_or_sce$end)
    )

  chr <- unique(ad_or_sce$chr)

  if(length(chr) != 1) {
    stop("The Credible sets are located on differenct chromosomes.")
  }

  idx <- which(
        rowData(ad_or_sce)$chr==chr & 
        rowData(ad_or_sce)$pos>=boundries["start"] & 
        rowData(ad_or_sce)$pos<=boundries["end"]
        )

  ad_or_sce <- ad_or_sce[idx,]

  return(ad_or_sce)

}