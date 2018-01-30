summariseArray <- function(mat, ids, method="rowmax"){
  
  dimnames(mat) -> dm
  
  if(method=="rowmax") by(mat, ids, function(x)x[which.max(rowMeans(x)),]) -> mat.by
  if(method=="average") by(mat, ids, function(x)colMeans(x)) -> mat.by
  t(sapply(mat.by, function(x)x)) -> mat.u
  apply(mat.u, 2, unlist) -> mat.u
  
  return(mat.u)
  
}