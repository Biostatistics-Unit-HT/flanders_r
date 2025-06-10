extract.abf.list.matrix=function(data=finemap.data,clique.list=structure_level){

  if(length(clique.list)>1){
    final.groups=clique.list[[length(clique.list)]]$clique.map
    for( livello in length(clique.list):2){
      for( i in 1:length(final.groups)){

        final.groups[[i]]=unique(unlist(clique.list[[livello-1]]$clique.map[final.groups[[i]]]))

      }
    }
    clique.list=final.groups


  }else{
    clique.list=clique.list[[1]]$clique.map
  }


  abf.clique=list()
  tmp.matr=t(data)
  abf.val=as.matrix(tmp.matr)

  for(mod in 1:length(clique.list)){

    if(length(clique.list[[mod]])>1){

      abf.clique[[mod]]=colSums(abf.val[clique.list[[mod]],])
      names(abf.clique[[mod]])=row.names(data)
    }else{
      abf.clique[[mod]]=abf.val[clique.list[[mod]],]
      names(abf.clique[[mod]])=row.names(data)

    }
  }

  return(abf.clique)

}
