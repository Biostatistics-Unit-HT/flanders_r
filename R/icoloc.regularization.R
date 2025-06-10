logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

logbf_to_pp=function(bf=all.coloc.join$ABF.tot) {
  denom=logsum(bf)
  exp(bf  - denom)
}



icoloc.regularization=function(cs.vector){
  cs.pip=logbf_to_pp(cs.vector)
  cs.vector=cbind(1:length(cs.vector),cs.vector,cs.pip,FALSE)
  cs.vector=cs.vector[order(cs.vector[,3],decreasing=T),]
  last.idx=min(which(cumsum(cs.vector[,3])>0.99))
  cs.vector[1:last.idx,4]=1
  cs.vector[cs.vector[,4]==0,2]=min(cs.vector[cs.vector[,4]==0,2],na.rm=T)
  cs.vector=cs.vector[order(cs.vector[,1],decreasing=F),]
  cs.vector[,2]
}
