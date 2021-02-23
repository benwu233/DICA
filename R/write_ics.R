#' @title Write the estimated components into a NIfTI file.
#'
#' @param ICs The output of DICA. A
#' @param nii The reference image.
#' @param mask A 3-D mask.
#' @param q The quantile thresholding the components.
#' @file Where the NIfTI file is saved to.
#'
#' @importFrom stats cov
#' @importFrom stats quantile
#' @importFrom neurobase copyNIfTIHeader
#' @importFrom oro.nifti writeNIfTI
#' @return
#' @export
#'
#' @examples
write_ics = function(ICs,nii,mask, q = 0.99, file){

  out = matrix(0, ncol = ncol(ICs), nrow = nrow(ICs))

  for(i in 1:ncol(ICs)){
    absIC = abs(ICs[,i])
    out[,i] = ICs[,i] * (absIC > quantile(absIC,q) )
  }

  dim0 = dim(mask)
  out_S = array(NaN, dim = c(dim0[1], dim0[2], dim0[3], ncol(ICs)))
  xgrid = as.matrix(expand.grid(1:dim0[1],1:dim0[2],1:dim0[3]))

  tag = 1

  for(i in 1:nrow(xgrid)){
    if(mask[xgrid[i,1],xgrid[i,2],xgrid[i,3]]==1){
      out_S[xgrid[i,1],xgrid[i,2],xgrid[i,3],] = out[tag,]
      tag = tag + 1
    }
  }

  copyNIfTIHeader(img = nii, arr = out_S)

  writeNIfTI(out_S,filename = file)
}
