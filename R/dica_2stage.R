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
#' @param Rmat Initial estimate of the \code{L*L} orthogonal rotation matrix at stage two.
#' @param tol Relative convergence tolerance for the log-likelihood.
#' @param itr.max Maximum number of algorithm iterations.
#' @param log if \code{TRUE}, the posterior probabilities at stage one are calculated on the log scale.
#' @param verbose If \code{TRUE}, print progress of algorithm to console.
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
#' @importFrom stats dnorm
#'
#' @examples
dica = function(Y, K=6, L=3, Rmat = diag(L),tol=1e-5,itr.max = 1000,log=FALSE,verbose = TRUE){

  cl0 = class(Y)[1]

  if(cl0 == "matrix"){
    J = nrow(Y)
    d = ncol(Y)

    if(verbose){
      print("Stage 1: Fitting a Mixture Distribution")
    }

    model1 = Mclust(Y, G = K, modelNames = "EEI",
                    control = emControl(tol = tol),verbose=verbose)
    #pi_j0 = model1$z

    logpi_j = matrix(0,nrow = J,ncol = K)
    pi_j = matrix(0,nrow = J,ncol = K)

    comp_sd = sqrt(diag(model1$parameters$variance$Sigma))
    comp_mean = model1$parameters$mean
    for(k in 1:K){
      tmp1 = dnorm(t(Y), mean = comp_mean[,k], sd = comp_sd,log=TRUE)
      logpi_j[,k] = apply(tmp1,2,sum)
    }

    if(!log){
      for(j in 1:J){
        tmp = exp(logpi_j[j,] - max(logpi_j[j,]))
        sumtmp = sum(tmp)
        pi_j[j,] = tmp / sumtmp
      }
    }

    order0 = order(apply(model1$parameter$mean,2,mean) )
    logpi_j = logpi_j[,order0]

    if(!log){
      pi_j = pi_j[,order0]
    }

    out = list()

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
    logpi_j = res$logposter[,order0]

    out = list()

  }
  else{
    return( print("The input is either a matrix or a list!") )
  }

  if(!log){
    mlogitProb = Mlogit(pi_j,exp(-300))
    if( sum(logpi_j < (-300) ) > 0.01*K*nrow(logpi_j) ){
      print("Warning: underflows occur. Consider computing on log scale with log=TRUE")
    }
  }else{
    mlogitProb = Mlogit2(logpi_j)
  }

  if(verbose){
    print("Stage 2: ICA Decomposition.")
  }
  res1 = icaimax(mlogitProb,L,center=FALSE,maxit = itr.max, Rmat = Rmat )
  S = res1$S
  A = res1$M

  out$S = S
  out$A = A

  if(verbose){
    print("Done.")
  }

  return(out)
}

