# compute graph with target number of edges

compute_graph_ziln <- function(X, target_edge_number, tol = 1) {
  l_min = 0
  l = 0.5
  l_max = 1
  res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
  edge_number = sum(res$path[[1]] != 0) /2
  while(edge_number > target_edge_number + tol || edge_number < target_edge_number - tol){
    if(edge_number > target_edge_number + tol) {
      l_ = (l_max + l) / 2
      l_min = l
      l = l_
    }
    else{
      l_ = (l_min + l) / 2
      l_max = l
      l = l_
    }
    print(l)
    res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
    edge_number = sum(res$path[[1]] != 0) /2
  }
  return(res$path[[1]])
}
