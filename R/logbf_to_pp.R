logbf_to_pp <- function(bf=all.coloc.join$ABF.tot) {
  denom=coloc:::logsum(bf)
  exp(bf  - denom)
}
