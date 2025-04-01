#' Splines for two-phase regression functions
#'
#' Creates splines for two-phase regression function in this package, including \code{linear2ph}, \code{logistic2ph}, \code{cv_linear2ph}, \code{cv_logistic2ph}.
#'
#' @param x  Column names of the covariate of the dataset.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param size Pass on to the \code{df} argument in \code{splines::bs()}. Degrees of freedom for EACH variable.
#' @param degree Pass on to the \code{degree} argument in \code{splines::bs()}. Degree of the piecewise polynomial. Default is 3 for cubic splines.
#' @param bs_names Column names of the output B-spline basis matrix.
#' @param group Optional. Column name of the categorical variable of which might have heterogeneous errors among different groups.
#' @param split_group Optional. Whether to split by group proportion for the group with B-spline size if the \code{group} argument is provided. If `FALSE`, then the split will be averaged across all groups. Default is `TRUE`.
#'
#' @details
#' This function can be directly applied for regression model with one or more error-prone continuous covariates.
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
spline2ph <- function(x, data, size = 20, degree = 3, bs_names, group = NULL, split_group = TRUE){
  n <- nrow(data) # get the number of rows of the dataset
  # portion the size of the B-spline basis according to the size of two sex groups
  if(is.null(group)){
    sn_all <- size
  } else {
    if(split_group){
      var_ratio <- table(data[[group]])/n
      group.names <- names(var_ratio)

    } else {
      group.names <- unique(data[[group]])
      n.group <- length(group.names)
      var_ratio <- rep((1/n.group), n.group)

    }
    sn_all <- floor(size*as.numeric(var_ratio))
    if(min(sn_all) < degree){stop("Size of B-spline is too small!")}
    if(size>sum(sn_all)){sn_all[1:(size-sum(sn_all))] <- sn_all[1:(size-sum(sn_all))]+1}
  }
  Bspline <- list()
  for(i in 1:length(sn_all)){
    # create B-spline basis for each variable group separately
    if(!is.null(group)){
      group.idx <- which(data[[group]]==group.names[i])
    }else{group.idx <- 1:n}

    sn_all.i <- sn_all[i] # size for the variable group
    BSpline.i <- matrix(data=0, nrow=n, ncol=(sn_all.i)^length(x)) # output matrix
    # create B-spline basis for each variable separately
    # then use tensor product to combine them
    idx.splines <- expand.grid(rep(list(1:sn_all.i), length(x)))
    names(idx.splines) <- x
    BSpline.x <- list()
    for(x.i in x){
      BSpline.x[[x.i]] <- splines::bs(x=data[[x.i]][group.idx], sn_all.i,
                                      degree=degree, intercept=TRUE)
    }
    for(bs.i in 1:nrow(idx.splines)){
      bs.val <- rep(1, length(group.idx))
      for(x.i in x){
        bs.val <- bs.val * BSpline.x[[x.i]][, idx.splines[bs.i, x.i]]
      }
      BSpline.i[group.idx, bs.i] <- bs.val
    }

    Bspline[[i]] <- BSpline.i
  }
  Bspline.df <- do.call(cbind, Bspline)


  # create B-spline basis for this analysis by combining the bases created
  Bspline.df <- data.frame(Bspline.df)
  colnames(Bspline.df) <- bs_names
  Bspline.df <- cbind(data, Bspline.df)

  return(Bspline.df)
}
