coloc.abf.matrix=function(clique.abfs=list(),threshold=0.75){
  W_comm=matrix(nrow = length(clique.abfs),ncol=length(clique.abfs))
  for(i in 1:(nrow(W_comm)-1)){
    for(j in (i+1):nrow(W_comm)){

      res=icoloc.ht(
        df1=data.frame(snp=names(clique.abfs[[i]]),lABF=unlist(clique.abfs[[i]])),
        df2=data.frame(snp=names(clique.abfs[[j]]),lABF=unlist(clique.abfs[[j]])),
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = 1e-5
      )
      W_comm[i,j]<-W_comm[j,i]<-res$summary$PP.H4.abf


    }

  }
  adiacency.matrix=Matrix::tril(W_comm)
  diag(adiacency.matrix)=0
  adiacency.matrix[adiacency.matrix<threshold]=0
  adiacency.matrix[adiacency.matrix>0]=1

  res=list(coloc.matrix=W_comm,adiacency.matrix=adiacency.matrix)
  return(res)
}
