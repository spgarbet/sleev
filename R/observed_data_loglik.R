# Helper function to convert column names/indices to numeric indices
get_col_indices <- function(cols, colnames_data)
{
  if (is.numeric(cols))   as.integer(cols)                       else
  if (is.character(cols)) as.integer(match(cols, colnames_data)) else
  stop("Columns must be numeric or character")
}

#' Observed-data log-likelihood
#'
#' This function returns the value of the observed-data log-likelihood (equation
#' (2) in Lotspeich et al. (2020)) for a given dataset and parameter values
#' `theta`, `gamma`, and `p`.
#'
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Predictors in the analysis model
#' @param gamma_pred Predictors in the outcome error model
#' @param theta Parameters for the analysis model (a column vector)
#' @param gamma Parameters for the outcome error model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function
#' @noRd

observed_data_loglik <- function(
  N,
  n,
  Y_unval      = NULL,
  Y            = NULL,
  X_unval      = NULL,
  X            = NULL,
  Z            = NULL,
  Bspline      = NULL,
  comp_dat_all,
  theta_pred,
  gamma_pred   = NULL,
  theta,
  gamma        = NULL,
  p)
{

  # Error settings based on arguments provided ---------------------------------
  errorsY <- !is.null(Y_unval) ## errors in Y
  errorsX <- !is.null(X_unval) ## errors in X

  m <- nrow(p)
  if(!is.matrix(comp_dat_all)) comp_dat_all <- data.matrix(comp_dat_all)
  colnames_data   <- colnames(comp_dat_all)
  theta_pred_cols <- get_col_indices(theta_pred, colnames_data)
  Y_col           <- get_col_indices(Y,          colnames_data)
  gamma_pred_cols <- get_col_indices(gamma_pred, colnames_data)
  Y_unval         <- get_col_indices(Y_unval,    colnames_data)
  bspline_cols    <- get_col_indices(Bspline,    colnames_data)
  k_col           <- get_col_indices("k",        colnames_data)

  ob_loglik(
    comp_dat_all    = comp_dat_all,
    N               = as.integer(N),
    n               = as.integer(n),
    theta_pred_cols = theta_pred_cols,
    theta           = as.numeric(theta),
    Y_col           = Y_col,
    gamma_pred_cols = gamma_pred_cols,
    gamma           = as.numeric(gamma),
    Y_unval_col     = Y_unval,
    bspline_cols    = bspline_cols,
    p               = p,
    k_col           = k_col
  )
}
