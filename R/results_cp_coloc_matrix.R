results_cp_coloc_matrix <- function(         W=Wmat
                                             ,      cliques=cliques_unweighted
                                             ,      labels=labels
                                             ,      finemap.data =finemap.data
){

  #retrieve abfs
  structure_level=list()
  structure_level[[1]]=list(clique.map=NA)
  structure_level[[1]]$clique.map=cliques
  level=1



  tmp.abf=extract.abf.list.matrix(data=finemap.data,clique.list=structure_level)
  print(6)

  W_comm=coloc.abf.matrix(clique.abfs = tmp.abf)
  W_comm=W_comm$adiacency.matrix
  fully.con=unweighted(qgraph::getWmat(qgraph::qgraph(W_comm)))
  W_i_comm <- igraph::graph_from_adjacency_matrix(as.matrix(W_comm),mode = "max")
  members <- igraph::components(W_i_comm)$membership
  print(7)



  while(length(fully.con)!=length(unique(members))){
    print(level)
    level=level+1
    structure_level[[level]]=list(clique.map=NA)
    structure_level[[level]]$clique.map=fully.con
    tmp.abf=extract.abf.list.matrix(data=finemap.data
                                    ,clique.list=structure_level)

    W_comm=coloc.abf.matrix(tmp.abf)
    W_comm=W_comm$adiacency.matrix
    fully.con=unweighted(qgraph::getWmat(qgraph::qgraph(W_comm)))
    W_i_comm <- igraph::graph_from_adjacency_matrix(as.matrix(W_comm))
    members <- (igraph::components(W_i_comm)$membership)
    if(length(fully.con)==length(unique(members))){
      level=level+1
      structure_level[[level]]=list(clique.map=NA)
      structure_level[[level]]$clique.map=fully.con

    }


  }

  final.groups=last(structure_level)$clique.map
  if(level>1){
    for( livello in length(structure_level):2){
      for( i in 1:length(final.groups)){

        final.groups[[i]]=unique(unlist(structure_level[[livello-1]]$clique.map[final.groups[[i]]]))

      }
    }
  }else{

    split=c()
    for(i in 1:length(fully.con)){

      split[[i]]=cliques[fully.con[[i]]]

    }

    #members <- igraph::components(W_i_comm)$membership
    #split <- split(cliques_unweighted, fully.con)
    final.groups <- lapply(split, function(x) sort(unique(unlist(x))))


  }


  print("adiacency created")

  communities=final.groups






  #if there is more than one community...
  #get shared nodes
  #if there are no shared nodes, create empty variable

  if (length(communities) > 1) {
    shared <- list()
    count <- 1
    for (i in 1:(length(communities) - 1)) {
      for (j in (i + 1):length(communities)) {
        shared[[count]] <- Reduce(intersect, list(communities[[i]],communities[[j]]))
        count <- count + 1
      }
    }
    shared <- sort(unique(unlist(shared)))
    if (length(shared) == 0) {shared <- c()}
  }

  #if there are communities...
  #create list of communities with labels from qgraph object instead of node numbers
  ### Add fix here

  if (length(communities) > 0) {
    communities_labels <- vector("list", length(communities))
    for (i in 1:length(communities)) {
      communities_labels[[i]] <- labels[communities[[i]]]
    }
  }

  #if there are no communities...
  if (length(communities) == 0) {
    communities_labels <- list()
  }


  #if there is more than one community...
  #create vector with shared nodes with labels from qgraph object instead of node numbers
  if (length(communities) > 1) {
    shared_labels <- labels[shared]
  }
  #if there is less than two communities
  if (length(communities) < 2) {
    shared <- c()
    shared_labels <- c()
  }
  #if there are no shared nodes, create empty variable
  if (length(shared_labels) == 0) {shared_labels <- c()}

  #return community lists, shared nodes vectors, and isolated nodes vectors
  return(list(communities,communities_labels,
              shared,shared_labels,
              W_comm))
}
