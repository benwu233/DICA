#' @title Reconstrct DTI observations.
#'
#' @param V1
#' @param V2
#' @param V3
#' @param L1
#' @param L2
#' @param L3
#'
#' @return
#' @export
#'
#' @examples
get_dti = function(V1,V2,V3,L1,L2,L3){
  dim0 = dim(V1)
  xgrid = as.matrix(expand.grid(1:dim0[1],1:dim0[2],1:dim0[3]))

  X = list()
  xgrid0 = NULL
  tag = 1
  for(i in 1:nrow(xgrid)){
    l1 = xgrid[i,1]
    l2 = xgrid[i,2]
    l3 = xgrid[i,3]
    if(L1[l1,l2,l3]!=0){
      xgrid0 = rbind(xgrid0,xgrid[i,])
      lambda = diag(c(L1[l1,l2,l3],L2[l1,l2,l3],L3[l1,l2,l3]))
      if(lambda[3,3] < 0){
        diag(lambda) = diag(lambda) - lambda[3,3]+ 1e-10
      }
      U = cbind(V1[l1,l2,l3,],V2[l1,l2,l3,],V3[l1,l2,l3,])
      tmp = U %*% lambda %*% t(U)
      X[[tag]] = (tmp + t(tmp) )/2
      tag = tag + 1
    }
  }

  out = list()

  out$X=X
  out$grid = xgrid0

  return(out)
}
