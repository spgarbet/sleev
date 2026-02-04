Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List e_step_complete_cpp(NumericMatrix theta_design_mat, NumericVector theta,
                         NumericMatrix comp_dat_all, int Y_col,
                         NumericMatrix gamma_design_mat, NumericVector prev_gamma,
                         int Y_unval_col, NumericMatrix prev_p,
                         IntegerVector Bspline_cols, int n, int N, int m,
                         bool errorsY)
{
  // Get the .pYstarCalc function from R
  Function pYstarCalc(".pYstarCalc");

  int N_minus_n = N - n;

  // Calculate pY_X using existing .pYstarCalc
  NumericVector pY_X = pYstarCalc(theta_design_mat, 1, theta,
                                   comp_dat_all, Y_col);

  // Calculate pYstar: P(Y*|X*,Y,X)
  NumericVector pYstar;
  if (errorsY)
  {
    pYstar = pYstarCalc(gamma_design_mat, n + 1, prev_gamma,
                       comp_dat_all, Y_unval_col);
  }
  else
  {
    pYstar = NumericVector(pY_X.size(), 1.0);
  }

  // Calculate pX: P(X|X*)
  // Extract B-spline columns
  int n_bspline_cols = Bspline_cols.size();
  NumericMatrix bspline_data(comp_dat_all.nrow() - n, n_bspline_cols);
  for (int i = 0; i < bspline_data.nrow(); ++i)
  {
    for (int j = 0; j < n_bspline_cols; ++j)
    {
      bspline_data(i, j) = comp_dat_all(n + i, Bspline_cols[j]);
    }
  }

  // Create pX matrix
  int pX_rows;
  if (errorsY)
  {
    pX_rows = m * N_minus_n * 2;
  }
  else
  {
    pX_rows = m * N_minus_n;
  }

  NumericMatrix pX(pX_rows, n_bspline_cols);

  // Fill pX with reordered prev_p multiplied by B-spline data
  if (errorsY)
  {
    // rep(rep(seq(1, m), each = (N - n)), times = 2)
    for (int rep_idx = 0; rep_idx < 2; ++rep_idx)
    {
      for (int k = 0; k < m; ++k)
      {
        for (int i = 0; i < N_minus_n; ++i)
        {
          int pX_row = rep_idx * m * N_minus_n + k * N_minus_n + i;
          for (int j = 0; j < n_bspline_cols; ++j)
          {
            pX(pX_row, j) = prev_p(k, j) * bspline_data(i, j);
          }
        }
      }
    }
  }
  else
  {
    // rep(seq(1, m), each = (N - n))
    for (int k = 0; k < m; ++k)
    {
      for (int i = 0; i < N_minus_n; ++i)
      {
        int pX_row = k * N_minus_n + i;
        for (int j = 0; j < n_bspline_cols; ++j)
        {
          pX(pX_row, j) = prev_p(k, j) * bspline_data(i, j);
        }
      }
    }
  }

  // Now calculate conditional expectations
  int n_rows = pX.nrow();
  int n_cols = pX.ncol();

  // Allocate matrices
  NumericMatrix psi_num(n_rows, n_cols);
  NumericMatrix psi_t(n_rows, n_cols);
  NumericVector w_t(n_rows);
  NumericMatrix u_t(m * N_minus_n, n_cols);

  // Calculate psi_num: c(pY_X * pYstar) * pX
  int pYstar_len = pYstar.size();

  for (int i = 0; i < n_rows; ++i)
  {
    // Handle R-style recycling for pYstar
    int pYstar_idx = (pYstar_len == 1) ? 0 : i;
    double weight = pY_X[i] * pYstar[pYstar_idx];

    for (int j = 0; j < n_cols; ++j)
    {
      psi_num(i, j) = weight * pX(i, j);
    }
  }

  // Calculate psi_denom by summing over groups using rowsum logic
  NumericVector psi_denom(N_minus_n, 0.0);

  for (int row = 0; row < n_rows; ++row)
  {
    int group_idx = row % N_minus_n;

    for (int j = 0; j < n_cols; ++j)
    {
      psi_denom[group_idx] += psi_num(row, j);
    }
  }

  // Avoid division by zero
  for (int i = 0; i < N_minus_n; ++i)
  {
    if (psi_denom[i] == 0.0)
    {
      psi_denom[i] = 1.0;
    }
  }

  // Calculate psi_t: psi_num / psi_denom
  for (int row = 0; row < n_rows; ++row)
  {
    int group_idx = row % N_minus_n;
    for (int j = 0; j < n_cols; ++j)
    {
      psi_t(row, j) = psi_num(row, j) / psi_denom[group_idx];
    }
  }

  // Calculate w_t by summing across columns (rowSums)
  for (int i = 0; i < n_rows; ++i)
  {
    w_t[i] = 0.0;
    for (int j = 0; j < n_cols; ++j)
    {
      w_t[i] += psi_t(i, j);
    }
  }

  // Calculate u_t
  int u_rows = m * N_minus_n;

  for (int i = 0; i < u_rows; ++i)
  {
    for (int j = 0; j < n_cols; ++j)
    {
      u_t(i, j) = psi_t(i, j);
      if (errorsY)
      {
        u_t(i, j) += psi_t(i + u_rows, j);
      }
    }
  }

  return List::create(
    Named("psi_t") = psi_t,
    Named("w_t") = w_t,
    Named("u_t") = u_t
  );
}
')



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
#' @importFrom stats as.formula
#' @importFrom stats glm
#'
#' @noRd

profile_out <- function(theta, n, N, Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL, Bspline = NULL,
                           comp_dat_all, theta_pred, gamma_pred, gamma0, p0, p_val_num, TOL, MAX_ITER)
{
  # Determine error settings
  errorsX <- !is.null(X_unval)
  errorsY <- !is.null(Y_unval)

  sn <- ncol(p0)
  m <- nrow(p0)

  prev_gamma <- gamma0
  prev_p <- p0

  # Convert to matrices
  theta_design_mat <- as.matrix(cbind(int = 1,
                                      comp_dat_all[-c(1:n), theta_pred]))
  comp_dat_all <- as.matrix(comp_dat_all)

  # Prepare column indices (convert to 0-based for C++)
  Y_col_idx <- match(Y, colnames(comp_dat_all)) - 1
  Y_unval_col_idx <- if (errorsY) match(Y_unval, colnames(comp_dat_all)) - 1 else 0L
  Bspline_col_idx <- match(Bspline, colnames(comp_dat_all)) - 1

  if (errorsY)
  {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- as.matrix(cbind(int = 1,
                                        comp_dat_all[, gamma_pred]))
  }
  else
  {
    gamma_design_mat <- matrix(0, 0, 0)  # Empty matrix for C++
  }

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1

  # pre-allocate memory for M-step variables
  if (errorsY)
  {
    mus_gamma <- vector("numeric", nrow(gamma_design_mat) * ncol(prev_gamma))
  }

  # Estimate gamma/p using EM
  while(it <= MAX_ITER && !CONVERGED)
  {
     ###################################################################
    # E Step
    e_step_results <- e_step_complete_cpp(
      theta_design_mat, theta, comp_dat_all, Y_col_idx,
      gamma_design_mat, prev_gamma, Y_unval_col_idx,
      prev_p, Bspline_col_idx, n, N, m, errorsY
    )

    psi_t <- e_step_results$psi_t
    w_t <- e_step_results$w_t
    u_t <- e_step_results$u_t

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    if (errorsY) {
      w_t <- .lengthenWT(w_t, n)
      muVector <- .calculateMu(gamma_design_mat, prev_gamma)
      gradient_gamma <- .calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[,c(Y_unval)], muVector)
      hessian_gamma <- .calculateHessian(gamma_design_mat, w_t, muVector, n, mus_gamma)

      new_gamma <- tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
                            error = function(err) {
                              matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (any(is.na(new_gamma))) {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula,
                                                 family = "binomial",
                                                 data = data.frame(comp_dat_all),
                                                 weights = w_t)$coefficients,
                                             ncol = 1))
      }

      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
    } else {
      new_gamma <- NULL
      gamma_conv <- TRUE
    }
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N ---------
    new_p_num <- p_val_num +
      rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p <- t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv <- abs(new_p - prev_p) < TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # End M Step
    ###################################################################

    all_conv <- c(gamma_conv, p_conv)
    if (mean(all_conv) == 1) { CONVERGED <- TRUE }

    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_gamma <- new_gamma
    prev_p <- new_p
    #  ------------------------------- Update values for next iteration

  }

  if(it == MAX_ITER & !CONVERGED)
  {
    CONVERGED_MSG <- "MAX_ITER reached"
    if (errorsY)
    {
      new_gamma <- matrix(data = NA,
                          nrow = nrow(gamma0),
                          ncol = 1)
    } else
    {
      new_gamma <- NULL
    }
    new_p <- matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
