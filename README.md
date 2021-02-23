
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DICA

<!-- badges: start -->
<!-- badges: end -->

## Installation

Install DICA with:

``` r
devtools::install_github("benwu233/DICA")
```

## Example

The package includes example data for fMRI and DTI analysis.

``` r
data_path = system.file("extdata",package="DICA") 
write_path = "/path/to/output/"
```

### fMRI

``` r
fmridata = oro.nifti::readNIfTI(paste0(data_path,"/fmri.nii.gz")) 
mask = 1 - is.na(fmridata[,,,1])

prep_data = prepfmri(fmridata,mask,40)
res_fmri = dica(prep_data,K = 20,L =14)

write_ics(res_fmri$S,fmridata,mask, q = 0.95,paste0(write_path,"fmri"))
```

### DTI

``` r
L_1 = oro.nifti::readNIfTI(paste0(data_path,"/dti_l1.nii.gz"))
L_2 = oro.nifti::readNIfTI(paste0(data_path,"/dti_l2.nii.gz"))
L_3 = oro.nifti::readNIfTI(paste0(data_path,"/dti_l3.nii.gz"))
V_1 = oro.nifti::readNIfTI(paste0(data_path,"/dti_v1.nii.gz"))
V_2 = oro.nifti::readNIfTI(paste0(data_path,"/dti_v2.nii.gz"))
V_3 = oro.nifti::readNIfTI(paste0(data_path,"/dti_v3.nii.gz"))

dtidata = get_dti(V_1,V_2,V_3,L_1,L_2,L_3)
res_dti = dica(dtidata$X,K = 20,L =14,tol = 1e-3)
write_ics(res_dti$S,L_1,(L_1!=0), q = 0.95, paste0(write_path,"dti"))
```