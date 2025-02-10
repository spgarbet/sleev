#' Splines for two-phase regression functions
#'
#' Supplements splines for two-phase regression function in this package, including \code{linear2ph}, \code{logistic2ph}, \code{cv_linear2ph}, \code{cv_logistic2ph}.
#'
#' @param x  Column name of the covariate of the dataset.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param size Pass on to the \code{df} argument in \code{splines::bs()}. Degrees of freedom.
#' @param degree Pass on to the \code{degree} argument in \code{splines::bs()}. Degree of the piecewise polynomial. Default is 3 for cubic splines.
#' @param group Optional. Column name of the categorical variable of which might have heterogeneous errors among different groups.
#'
#' @details
#' This function can be directly applied to regression model with one error-prone covariate. If there are more than one error-prone covariate, please check vignette on how to handle it.
#'
#' @return
#' the \code{data.frame} object including the original dataset and the B-spline bases.
#'
#' @examples
#' # example code
#' data("mock.vccc")
#' spline2ph(size = 20, degree = 3, data = mock.vccc, x = "VL_unval", group = "Sex")
#'
#' @export
spline2ph <- function(size = 20, degree = 3, data, x, group = NULL){
  n <- nrow(data) # get the number of rows of the dataset
  # portion the size of the B-spline basis according to the size of two sex groups
  var_ratio <- table(data[[group]])/n
  group.names <- names(var_ratio)
  sn_all <- round(size*as.numeric(var_ratio), digits = 0)

  Bspline <- list()
  for(i in 1:length(sn_all)){
    group.idx <- which(data[[group]]==group.names[i])
    # create B-spline basis for each sex group separately through bs()
    BSpline.i <- matrix(data=0, nrow=n, ncol=sn_all[i])
    BSpline.i[group.idx, ]<- splines::bs(x=data[[x]][group.idx], df=sn_all[i],
                degree=degree, intercept=TRUE)
    Bspline[[i]] <- BSpline.i
  }
  # create B-spline basis for this analysis by combining the bases created
  Bspline.df <- do.call(cbind, Bspline)
  colnames(Bspline.df) <- paste0("bs", 1:size)
  new.data <- cbind(data, Bspline.df)
  return(new.data)
}
