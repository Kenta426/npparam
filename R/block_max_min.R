block <- function(y){
  # y_{ij} = f(x_{ij}) + eps_{ij}, where {x_{ij}} form a lattice
  
  ## 1. pre-processing
  n <- prod(dim(y)) # total number of points
  V <- which(array(rep(1,n), dim(y))==1, arr.ind = T) # vector of vertices
  A <- t(sapply(1:n, FUN = function(k){
    return(sapply(1:n, FUN =function(l){
      return(prod(V[k,] <= V[l,]))}))})) # adjacency matrix
  n1 <- dim(y)[1]
  E.coor <- which(A==1, arr.ind = T)
  E <- cbind((E.coor - 1) %% n1 + 1, (E.coor - 1) %/% n1 + 1) # vector of edges, represented by (col1, col2) -> (col3, col4) 
  E_index <- which(A==1) # vector of edges, represented by the indices of nonzero adjacency matrix entries
  
  ## 2. calculate block average matrix mu, mu_[u,v] = Average of y_{ij} for {x_{ij} in [u,v]}
  mu = matrix(nrow = n, ncol = n)
  mu[E_index] = apply(E, MARGIN = 1, FUN = function(vec) { 
    return(mean(y[vec[1] : vec[2], vec[3]:vec[4]]))
  })
  
  ## 3. calculate the block max-min estimator
  output = lapply(1:n, FUN = function(j){
    u.ind = which(A[,j]==1)
    v.ind = which(A[j,]==1)
    min_output = sapply(u.ind, function(x){
      return(min(mu[x, v.ind]))
    })
    return(max(min_output))
  })
  
  ## 4. return value
  return(array(unlist(output), dim = dim(y)))
}
