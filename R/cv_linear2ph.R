#' Performs cross-validation to calculate the average predicted log likelihood for the \code{linear2ph} function. This function can be used to select the B-spline basis that yields the largest average predicted log likelihood.
#'
#' @param y_unval Specifies the column of the error-prone outcome that is continuous. Subjects with missing values of \code{Y_unval} are omitted from the analysis. This argument is required.
#' @param y Specifies the column that stores the validated value of \code{Y_unval} in the second phase. Subjects with missing values of \code{Y} are considered as those not selected in the second phase. This argument is required.
#' @param x_unval Specifies the columns of the error-prone covariates. Subjects with missing values of \code{X_unval} are omitted from the analysis. This argument is required.
#' @param x Specifies the columns that store the validated values of \code{X_unval} in the second phase. Subjects with missing values of \code{X} are considered as those not selected in the second phase. This argument is required.
#' @param b_spline Specifies the columns of the B-spline basis. Subjects with missing values of \code{Bspline} are omitted from the analysis. This argument is required.
#' @param z Specifies the columns of the accurately measured covariates. Subjects with missing values of \code{Z} are omitted from the analysis. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param nfolds Specifies the number of cross-validation folds. The default value is \code{5}. Although \code{nfolds} can be as large as the sample size (leave-one-out cross-validation), it is not recommended for large datasets. The smallest value allowable is \code{3}.
#' @param max_iter Specifies the maximum number of iterations in the EM algorithm. The default number is \code{2000}. This argument is optional.
#' @param tol Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#' @return
#' \item{avg_pred_loglike}{Stores the average predicted log likelihood.}
#' \item{pred_loglike}{Stores the predicted log likelihood in each fold.}
#' \item{converge}{Stores the convergence status of the EM algorithm in each run.}
#' @importFrom Rcpp evalCpp
#' @importFrom stats pchisq
#' @examples
#'   rho = 0.3
#'   p = 0.3
#'   n = 100
#'   n2 = 40
#'   alpha = 0.3
#'   beta = 0.4
#'
#'   ### generate data
#'   simX = rnorm(n)
#'   epsilon = rnorm(n)
#'   simY = alpha+beta*simX+epsilon
#'   error = MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1, rho, rho, 1), nrow=2))
#'
#'   simS = rbinom(n, 1, p)
#'   simU = simS*error[,2]
#'   simW = simS*error[,1]
#'   simY_tilde = simY+simW
#'   simX_tilde = simX+simU
#'
#'   id_phase2 = sample(n, n2)
#'
#'   simY[-id_phase2] = NA
#'   simX[-id_phase2] = NA
#'
#'   # cubic basis
#'   nsieves = c(5, 10)
#'   pred_loglike = rep(NA, length(nsieves))
#'   for (i in 1:length(nsieves)) {
#'       nsieve = nsieves[i]
#'       Bspline = splines::bs(simX_tilde, df=nsieve, degree=3,
#'         Boundary.knots=range(simX_tilde), intercept=TRUE)
#'       colnames(Bspline) = paste("bs", 1:nsieve, sep="")
#'       # cubic basis
#'
#'       data = data.frame(Y_tilde=simY_tilde, X_tilde=simX_tilde, Y=simY, X=simX, Bspline)
#'       ### generate data
#'
#'       res = cv_linear2ph(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde",
#'         Bspline=colnames(Bspline), data=data, nfolds = 5)
#'       pred_loglike[i] = res$avg_pred_loglik
#'     }
#'
#'   data.frame(nsieves, pred_loglike)
#'
#' @export
cv_linear2ph <- function (y_unval=NULL, y=NULL, x_unval=NULL, x=NULL, z=NULL, b_spline=NULL, data=NULL, nfolds=5, max_iter=2000, tol=1E-4, verbose=FALSE) {

  # variable name change
  Y_unval = y_unval; Y = y ; X_unval = x_unval; X = x; Z = z
  Bspline = b_spline; TOL = tol; MAX_ITER = max_iter

  ###############################################################################################################
  #### check data ###############################################################################################
  storage.mode(MAX_ITER) = "integer"
  storage.mode(TOL) = "double"
  storage.mode(nfolds) = "integer"

  if (missing(data)) {
    stop("No dataset is provided!")
  }

  if (missing(Y_unval)) {
    stop("The error-prone response Y_unval is not specified!")
  } else {
    vars_ph1 = Y_unval
  }

  if (missing(X_unval)) {
    stop("The error-prone covariates X_unval is not specified!")
  } else {
    vars_ph1 = c(vars_ph1, X_unval)
  }

  if (missing(Bspline)) {
    stop("The B-spline basis is not specified!")
  } else {
    vars_ph1 = c(vars_ph1, Bspline)
  }

  if (missing(Y)) {
    stop("The accurately measured response Y is not specified!")
  }

  if (missing(X)) {
    stop("The validated covariates in the second-phase are not specified!")
  }

  if (length(X_unval) != length(X)) {
    stop("The number of columns in X_unval and X is different!")
  }

  if (!missing(Z)) {
    vars_ph1 = c(vars_ph1, Z)
  }

  id_exclude = c()
  for (var in vars_ph1) {
    id_exclude = union(id_exclude, which(is.na(data[,var])))
  }

  if (verbose) {
    print(paste("There are", nrow(data), "observations in the dataset."))
    print(paste(length(id_exclude), "observations are excluded due to missing Y_unval, X_unval, or Z."))
  }
  if (length(id_exclude) > 0) {
    data = data[-id_exclude,]
  }

  n = nrow(data)
  if (verbose) {
    print(paste("There are", n, "observations in the analysis."))
  }

  id_phase1 = which(is.na(data[,Y]))
  for (var in X) {
    id_phase1 = union(id_phase1, which(is.na(data[,var])))
  }
  if (verbose) {
    print(paste("There are", n-length(id_phase1), "observations validated in the second phase."))
  }

  if (nfolds >= 3) {
    if (verbose) {
      print(paste0(nfolds, "-folds cross-validation will be performed."))
    }
  } else {
    stop("nfolds needs to be greater than or equal to 3!")
  }
  #### check data ###############################################################################################
  ###############################################################################################################



  ###############################################################################################################
  #### prepare analysis #########################################################################################
  Y_unval_vec = c(as.vector(data[-id_phase1,Y_unval]), as.vector(data[id_phase1,Y_unval]))
  storage.mode(Y_unval_vec) = "double"

  X_unval_mat = rbind(as.matrix(data[-id_phase1,X_unval]), as.matrix(data[id_phase1,X_unval]))
  storage.mode(X_unval_mat) = "double"

  Bspline_mat = rbind(as.matrix(data[-id_phase1,Bspline]), as.matrix(data[id_phase1,Bspline]))
  storage.mode(Bspline_mat) = "double"

  Y_vec = as.vector(data[-id_phase1,Y])
  storage.mode(Y_vec) = "double"

  X_mat = as.matrix(data[-id_phase1,X])
  storage.mode(X_mat) = "double"

  if (!is.null(Z)) {
    Z_mat = rbind(as.matrix(data[-id_phase1,Z]), as.matrix(data[id_phase1,Z]))
    storage.mode(Z_mat) = "double"
  }

  if (is.null(Z)) {
    Z_mat = rep(1., n)
  } else {
    Z_mat = cbind(1, Z_mat)
  }

  idx_fold = c(sample(1:nfolds, size = length(Y_vec), replace = TRUE),
               sample(1:nfolds, size = length(id_phase1), replace = TRUE))
  #### prepare analysis #########################################################################################
  ###############################################################################################################



  ###############################################################################################################
  #### analysis #################################################################################################
  pred_loglik = rep(NA, nfolds)
  converge = rep(NA, nfolds)
  for (fold in 1:nfolds) {
    Train = as.numeric(idx_fold != fold)
    res = .TwoPhase_MLE0_MEXY_CV_loglik(Y_unval_vec, X_unval_mat, Y_vec, X_mat, Z_mat, Bspline_mat, MAX_ITER, TOL, Train)
    pred_loglik[fold] = res$pred_loglike
    converge[fold] = !res$flag_nonconvergence
    if (pred_loglik[fold] == -999.) {
      pred_loglik[fold] = NA
    }
  }

  #### analysis #################################################################################################
  ###############################################################################################################



  ###############################################################################################################
  #### return results ###########################################################################################
  avg_pred_loglik = mean(pred_loglik, na.rm = TRUE)
  res_final = list(avg_pred_loglik=avg_pred_loglik, pred_loglik=pred_loglik, converge=converge)
  res_final
  #### return results ###########################################################################################
  ###############################################################################################################
}
