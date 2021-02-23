mixWishart = function(data,k=5,tol = 0.01,itr.max = 100){
  n = length(data)

  dataM = matrix(0, nrow = n, ncol = 9)

  for(i in 1:n){
    dataM[i,] = as.numeric(data[[i]])
  }

  tmp = kmeans(dataM,k,iter.max = 612)

  #initial
  V = array(0,dim=c(k,3,3))
  df = rep(5,k)

  for(j in 1:k){
    tmp1 = tmp$centers[j,]
    dim(tmp1) = c(3,3)
    V[j,,] = tmp1/df[j]
  }

  pi = rep(1/k,k)

  #tol0 = tol + 0.1
  tol0 = NULL

  tag = 0
  loglik = rep(0,itr.max)
  while(tag < itr.max){

    tag = tag  + 1

    #E-step

    pi_jk = matrix(0, nrow = n, ncol= k)
    pi_out = matrix(0, nrow = n, ncol= k)

    sumX = array(0, dim= c(3,3,k) )
    sumX_inv = array(0, dim= c(3,3,k) )
    sumlogX = rep(0,k)
    nk = rep(0,k)

    for(i in 1:n){

      for(j in 1:k){
        pi_jk[i,j] = dWISHART(data[[i]],df[j],V[j,,],log=TRUE) + log(pi[j])
      }

      ind = which.max(pi_jk[i,])
      sumX[,,ind] = sumX[,,ind] + data[[i]]
      nk[ind] = nk[ind] + 1
      sumlogX[ind] = sumlogX[ind]  + log( det(data[[i]])+1e-40 )

      loglik[tag] = loglik[tag] + pi_jk[i,ind]
    }

    for(j in 1:k){
      pi[j] = nk[j]/n
      sumX_inv[,,j] = solve(sumX[,,j] )
    }

    #M_step
    if(tag>1){
      tol0 = abs( (loglik[tag] - loglik[tag-1])/loglik[tag-1])
      if( tol0 < tol){

        for(i in 1:n){
          tmp = exp( pi_jk[i,] - max(pi_jk[i,]) )
          pi_out[i,] = tmp / sum(tmp)
        }

        loglik = loglik[1:tag]
        break
      }
    }

    V_old = V
    df_old = df
    for(j in 1:k){
      if(nk[j] > 0){
        tmp = mleWishart(nk[j],sumX_inv[,,j],sumlogX[j],solve(V_old[j,,]),(df_old[j]-4)/2,tol = 0.0001 )
        V[j,,] = tmp$V
        df[j] = tmp$df
      }
    }

    #print(paste0("loglik: ",loglik[tag]) )
    print(paste0("itr: ",tag, " tol: ", tol0))

  }

  out = list()
  out$V = V
  out$df = df
  out$pi = pi
  out$poster = pi_out
  out$loglik = loglik

  return(out)
}


mleWishart = function(N,sumX_inv,sumlogX,theta_s=NA,theta_df=NA,tol=0.01){

  tol0 = tol + 0.1
  tag = 0
  dim0 = nrow(sumX_inv)

  if(is.na(theta_s[1])){
    theta_s = diag(1,dim0)
  }

  if(is.na(theta_df)){
    theta_df = 1
  }

  while(tol0 > tol){

    tag = tag + 1
    theta_s_new = (2*theta_df + dim0 + 1) * N * sumX_inv

    tmp = sumlogX/N - dim0* log(2) + log( det(theta_s_new))

    gap = 100
    lower0 = 1+1e-10
    upper0 = 1e300
    if( as.numeric( mvdigamma(lower0,dim0) )> tmp){
      target = lower0
    }
    else if( as.numeric( mvdigamma(upper0,dim0) )< tmp){
      target = upper0
    }
    else{
      while(gap>0.00001){
        mid0 = (lower0 + upper0)/2
        tmp1 = as.numeric( mvdigamma(mid0,dim0) )

        if( tmp1 < tmp ){
          lower0 = mid0
        }
        else{
          upper0 = mid0
        }

        gap = abs(tmp1 - tmp)
      }

      target = mid0
    }

    theta_df_new = target - (dim0+1)/2

    tol0 = max(abs( (theta_df_new - theta_df)/(theta_df+1e-30) ), abs(theta_s_new - theta_s)/ (theta_s + 1e-30) )

    theta_df = theta_df_new
    theta_s = theta_s_new

    #print( paste0("itr: ", tag, " tol: ", tol0) )
  }

  out  = list()

  out$theta_s = theta_s
  out$theta_df = theta_df
  out$V = solve(theta_s)
  out$df = 2 * theta_df + dim0 + 1

  return(out)
}

Mlogit = function(Pi,e = 0.1){

  K = ncol(Pi)

  Pi = Pi+e
  out =  log( Pi[,-K] / Pi[,K]  )

  return(out)
}
