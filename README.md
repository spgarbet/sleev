<p style="display:inline-block;">
  <img src="hex.png" width="200">
    <h1>SLEEV: Semiparametric Likelihood Estimation with Errors in Variables</h1>
    </p>
    
The complete R package `logreg2ph` and code for the simulation settings included in [this paper](https://doi.org/10.1111/biom.13512). 
  
Part of the R package `TwoPhaseReg` as described in [this paper](https://doi.org/10.1002/sim.8876).
  
  This package combines elements of the R packages [`logreg2ph`](https://github.com/sarahlotspeich/logreg2ph) and [`TwoPhaseReg`](https://github.com/dragontaoran/TwoPhaseReg)
  
Below is a simple example for running a linear regression with `sleev`. The complete vignette detailing all uses of example in `sleev` can be found in folder `vignette`.
  
### Install
To install the package, run the following in your `R` console: 
    
```{r}
  devtools::install_github("dragontaoran/sleev")
```
  
Load package:
    
```{r}
  library(sleev)
```
  
### `mock.vccc` data
  
This data is a simulated two-phase sampling dataset based on VCCC data. To load this dataset, run:
    
```{r}
  data(mock.vccc)
```
  
### Preprocessing
  Because of skewness, we often transform both CD4 and VL. In our analysis, CD4 was divided by 10 and square-root transformed and VL was log_10 transformed:

```
  mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val/10)
  mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval/10)
  mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)
  mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)
```
Load `splines`	package
  
```
  library(splines)
```
  
### Construct B-spline Basis
  
```
  n <- nrow(mock.vccc) # get the number of rows of the dataset
  sn <- 20 # set the size of the B-spline basis
  # portion the size of the B-spline basis according to the size of two sex groups
  sex_ratio <- sum(mock.vccc$Sex==0)/n
  sn_0 <- round(sn*sex_ratio, digits = 0)
  sn_1 <- sn-sn_0
  # create B-spline basis for each sex group separately through bs()
  Bspline_0 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$Sex==0], df=sn_0, degree=3,
                           intercept=TRUE)
  Bspline_1 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$Sex==1], df=sn_1, degree=3,
                           intercept=TRUE)
  # create B-spline basis for this analysis by combining the two bases created above
  Bspline <- matrix(data=0, nrow=n, ncol=sn)
  Bspline[mock.vccc$Sex==0,1:(sn_0)] <- Bspline_0
  Bspline[mock.vccc$Sex==1,(sn_0+1):sn] <- Bspline_1
  # name B-spline basis matrix columns properly
  colnames(Bspline) <- paste0("bs", 1:sn) 
  # add the B-spline basis to the analysis dataset 
  data <- data.frame(cbind(mock.vccc, Bspline))
```
  
### Model fitting 
  
The SMLEs can be obtained by running
```
  res_linear <- linear2ph(Y_unval="CD4_unval_sq10", Y="CD4_val_sq10", X_unval="VL_unval_l10",
                          X="VL_val_l10", Z="Sex", Bspline=colnames(Bspline), data=data,
                          hn_scale = 1, noSE = FALSE, TOL = 1e-04, MAX_ITER = 1000,
                          verbose = FALSE)
```
Check convergence:
    
```
  c(res_linear$converge, res_linear$converge_cov)
```
  
