#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,Z), this function repeats the EM algorithm to find
#' the values of `gamma` and `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param gamma0 Starting values for `gamma`, the parameters for the outcome error model (a column vector)
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#'
#' @return Profile likelihood for `theta`: the value of the observed-data log-likelihood after profiling out other parameters.
#'
#' @noRd

profile_out <- function(theta, n, N, Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL, Bspline = NULL,
                        comp_dat_all, theta_pred, gamma_pred, gamma0, p0, p_val_num, TOL, MAX_ITER)
{
  # Determine error settings
  errorsX <- !is.null(X_unval)
  errorsY <- !is.null(Y_unval)

  m <- nrow(p0)

  # Convert to matrices
  comp_dat_all <- as.matrix(comp_dat_all)
  theta_design_mat <- as.matrix(cbind(int = 1,
                                      comp_dat_all[-c(1:n), theta_pred]))

  # Prepare column indices (convert to 0-based for C++)
  Y_col_idx <- match(Y, colnames(comp_dat_all)) - 1
  Y_unval_col_idx <- if (errorsY) match(Y_unval, colnames(comp_dat_all)) - 1 else 0L
  Bspline_col_idx <- match(Bspline, colnames(comp_dat_all)) - 1

  # Prepare gamma design matrix
  gamma_design_mat <- if (errorsY)
    as.matrix(cbind(int = 1, comp_dat_all[, gamma_pred])) else
    gamma_design_mat <- matrix(0, nrow = 0, ncol = 0) # Empty matrix

  # Call the complete EM loop in C++
  result <- profile_out_em(
    theta_design_mat = theta_design_mat,
    theta = theta,
    comp_dat_all = comp_dat_all,
    Y_col = Y_col_idx,
    gamma_design_mat = gamma_design_mat,
    gamma0 = gamma0,
    Y_unval_col = Y_unval_col_idx,
    p0 = p0,
    Bspline_cols = Bspline_col_idx,
    p_val_num = p_val_num,
    n = n,
    N = N,
    m = m,
    errorsY = errorsY,
    TOL = TOL,
    MAX_ITER = MAX_ITER
  )

  return(result)
}
