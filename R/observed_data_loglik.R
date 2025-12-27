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

#
#   ## Calculate P(Yi*|Xi*,y,xk) for all (y,xk) ------------------------------------
  if (errorsY)
  {
    pYstar <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), gamma_pred]) %*% gamma)))
    pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)]
   # cat("R pYstar", paste0(pYstar[1:5], collapse=", "), '\n')
  }
  ## ------------------------------------ Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  if (errorsX)
  {
    pX <- p[comp_dat_all[-c(1:n), "k"], ]
    #cat("R pX", paste0(pX[1,1:5], collapse=", "), '\n')
  }
#   ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
#
#   ################################################################################
#   # For unvalidated subjects ------------------------------------------------------
#   ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
    pY_X <- 1 / (1 + exp(-as.numeric((cbind(int = 1, comp_dat_all[-c(1:n), theta_pred]) %*% theta))))
    pY_X[which(comp_dat_all[-c(1:n), Y] == 0)] <- 1 - pY_X[which(comp_dat_all[-c(1:n), Y] == 0)]

    #cat("R pY_X", paste0(pY_X[1:5], collapse=", "), '\n')

    # if (errorsY & errorsX) {
    #   temp <- c(pY_X * pYstar * pX) * comp_dat_all[-c(1:n), Bspline]
    #   cat("temp dimensions:", dim(temp), "\n")
    #   cat("temp first 10:", head(c(temp), 10), "\n")
    #
    #   person_sum_mat <- rowsum(temp, group = rep(seq(1, (N - n)), times = 2 * m))
    #   cat("person_sum_mat dimensions:", dim(person_sum_mat), "\n")
    #   cat("First 3 rows of person_sum_mat:\n")
    #   print(person_sum_mat[1:3, ])
    #
    #   person_sum <- rowSums(person_sum_mat)
    #   cat("First 5 person_sum:", head(person_sum, 5), "\n")
    # }

#   ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
#   ################################################################################
#   ## Calculate sum of P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj --------------------
#   if (errorsY & errorsX)
#   {
#     person_sum <- rowsum(c(pY_X * pYstar * pX) * comp_dat_all[-c(1:n), Bspline],
#                          group = rep(seq(1, (N - n)),
#                                      times = 2 * m))
#   } else if (errorsY)
#   {
#     person_sum <- rowsum(c(pY_X * pYstar),
#                          group = rep(seq(1, (N - n)),
#                                      times = 2))
#   } else if (errorsX)
#   {
#     person_sum <- rowsum(c(pY_X * pX) * comp_dat_all[-c(1:n), Bspline],
#                          group = rep(seq(1, (N - n)),
#                                      times = m))
#   }
#   person_sum <- rowSums(person_sum)
#   log_person_sum <- log(person_sum)
#   log_person_sum[log_person_sum == -Inf] <- 0
#   ## And sum over them all -------------------------------------------------------
#
#   return_loglik + sum(log_person_sum)

    # After computing person_sum_mat
    # cat("\nDetailed person 0 computation:\n")
    # cat("person_sum_mat[1, ]:", person_sum_mat[1, ], "\n")
    # cat("sum(person_sum_mat[1, ]):", sum(person_sum_mat[1, ]), "\n")
    #
    # if (errorsY & errorsX) {
    #   # Show the grouping vector
    #   group_vec <- rep(seq(1, (N - n)), times = 2 * m)
    #   cat("Group vector length:", length(group_vec), "\n")
    #   cat("First 20 group values:", head(group_vec, 20), "\n")
    #   cat("Expected:", 2 * m, "rows per person\n")
    #
    #   # Create temp matrix
    #   temp <- c(pY_X * pYstar * pX) * comp_dat_all[-c(1:n), Bspline]
    #   cat("temp is a vector? ", is.vector(c(temp)), " length:", length(c(temp)), "\n")
    #
    #   cat("Rows for person 1:", which(group_vec == 1), "\n")
    #   cat("Rows for person 2:", which(group_vec == 2), "\n")
    # }

    ob_loglik(
      comp_dat_all    = comp_dat_all,
      N               = as.integer(N),
      n               = as.integer(n),
      theta_pred_cols = theta_pred_cols,
      theta           = as.numeric(theta),
      Y_col           = Y_col,
      gamma_pred_cols = gamma_pred_cols,
      gamma           = as.numeric(gamma),
      Y_unval         = Y_unval,
      bspline_cols    = bspline_cols,
      p               = p,
      k_col           = k_col
    )
}
