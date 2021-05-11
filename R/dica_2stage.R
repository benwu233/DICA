#' @title Distributional ICA
#' @description Distributional independent component ICA with two-stage estimation.
#' At stage one, a mixture model (Gaussian or Wishart) is estimated by EM algorithm. At stage two, the posterior
#' weights (after a transformation with the multinomial logit link function)
#' are decomposed with Infomax ICA.
#'
#' @param Y A numeric matrix or a list. If a matrix, a mixture of Gaussian model
#' is to be fitted at stage one, and rows correspond to observations (\code{n}) and columns
#' correspond to the dimension (\code{d}). If a list, a mixture of Wishart model
#' is to be fitted at stage one. Each item of the list is a \code{d*d} positive definite
#' matrix, corresponding to a observation.
#' @param K The number of components at stage one.
#' @param L The number of components at stage two, should be no larger than \code{K}.
#' @param tol Relative convergence tolerance for the log-likelihood.
#' @param itr.max Maximum number of algorithm iterations.
#' @param verbose If TRUE, print progress of algorithm to console.
#'
#' @return A list includes
#' \itemize{
#' \item \code{S} - The estimated component matrix.
#' \item \code{A} - The estiamted mixing matrix.
#' }
#' @export
#' @import mclust
#' @importFrom ica icaimax
#' @importFrom mixAK dWISHART
#' @importFrom CholWishart mvdigamma
#' @importFrom stats kmeans
#'
#' @examples
dica = function(Y, K=6, L=3,tol=1e-5,itr.max = 1000,verbose = TRUE){

  cl0 = class(Y)[1]

  if(cl0 == "matrix"){
    J = nrow(Y)
    d = ncol(Y)

    if(verbose){
      print("Stage 1: Fitting a Mixture Distribution")
    }

    model1 = Mclust(Y, G = K, modelNames = "EEI",
                    control = emControl(tol = tol),verbose=verbose)
    pi_j = model1$z

    order0 = order(apply(model1$parameter$mean,2,mean) )
    pi_j = pi_j[,order0]

  }
  else if(cl0 == "list"){

    if(verbose){
      print("Stage 1: Fitting a Mixture Distribution.")
    }

    res = mixWishart(Y,k = K,tol= tol,itr.max,verbose)

    det0 = rep(0,K)
    for(j in 1:K){
      det0[j] = det(res$V[j,,])
    }

    order0 = order(det0 )
    pi_j = res$poster[,order0]

  }
  else{
    return( print("The input is either a matrix or a list!") )
  }

  mlogitProb = Mlogit(pi_j,exp(-300))

  if(verbose){
    print("Stage 2: ICA Decomposition.")
  }
  res1 = icaimax(mlogitProb,L,center=FALSE,maxit = itr.max )
  S = res1$S
  A = res1$M

  out = list()

  out$S = S
  out$A = A

  if(verbose){
    print("Done.")
  }

  return(out)
}

