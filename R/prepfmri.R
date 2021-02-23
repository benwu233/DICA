#' @title preprocessing fMRI data
#' @description The function preprocesses fMRI data (NIfTI format) with the following steps:
#' detrending; demeaning; and PCA.
#'
#' @param nii A 4-D fMRI image.
#' @param mask A 3-D mask has the same dimension with \code{dim(nii)[1:3]}.
#' @param r The number of PCs.
#'
#'
#' @return The preprocessed n*r data matrix (n is the number of masked voxels).
#' @export
#'
#' @examples
prepfmri = function(nii,mask,r=40) {

  temp_data = matrix(0, nrow = dim(nii)[1]*dim(nii)[2]*dim(nii)[3], ncol = dim(nii)[4])
  temp_mask = rep(0, dim(mask)[1]*dim(mask)[2]*dim(mask)[3])

  xgrid = as.matrix(expand.grid(1:dim(nii)[1],1:dim(nii)[2],1:dim(nii)[3]))

  for(i in 1:nrow(xgrid)){
    temp_data[i,] = nii[xgrid[i,][1],xgrid[i,][2],xgrid[i,][3],]
    temp_mask[i] = mask[xgrid[i,][1],xgrid[i,][2],xgrid[i,][3]]
  }

  temp = temp_data[temp_mask==1,]

  J = nrow(temp)

  # detrend
  t = (1:dim(nii)[4]) - mean(1:dim(nii)[4])
  x = cbind(1,t,t^2,t^3)
  y = t(temp)
  beta = solve(t(x)%*% x) %*%t(x) %*% y
  y_hat = x%*%beta
  temp = t(y - y_hat)

  # demean
  temp_demeaned = matrix(0,nrow = nrow(temp), ncol = ncol(temp))
  meantemp = apply(temp,2,mean)
  for(i in 1:nrow(temp)) {
    temp_demeaned[i,] = temp[i,] - meantemp
  }

  # svd
  S_var = cov(temp_demeaned)
  egS_var = eigen(S_var)
  R = min(r,ncol(S_var) )

  U_R = egS_var$vectors[,1:R]
  lambda_R = diag(egS_var$values[1:R])

  reduced_data = temp_demeaned %*% U_R

  return(reduced_data)
}
