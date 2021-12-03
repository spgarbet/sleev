#' Cross-validated observed-data log-likelihood
#' This function returns the value of the observed-data log-likelihood based on cross-validation.
#'
#' @param fold Column name with the assigned fold for cross-validation.
#' @param Y_unval Column name with the unvalidated outcome. If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y Column name with the validated outcome.
#' @param X_unval Column name(s) with the unvalidated predictors.  If \code{X_unval} and \code{X} are \code{null}, all precictors are assumed to be error-free.
#' @param X Column name(s) with the validated predictors. If \code{X_unval} and \code{X} are \code{null}, all precictors are assumed to be error-free.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y}, \code{X_unval}, \code{X}, \code{Z}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return scalar value of the function
#' @export

cv_observed_data_loglik <- function(fold, Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL,
                        Validated = NULL, Bspline = NULL, data, theta_pred = NULL, gamma_pred = NULL,
                        TOL = 1E-4, MAX_ITER = 1000) {
  if (is.null(theta_pred)) { theta_pred <- c(X, Z) }
  if (is.null(gamma_pred) & !is.null(Y_unval)) { gamma_pred <- c(X_unval, Y, X, Z) }

  num_folds <- length(unique(data[, fold]))
  status <- rep(TRUE, num_folds)
  msg <- rep("", num_folds)
  ll <- rep(NA, num_folds)
  #fold_ll <- re_fold_ll <- vector()
  for (i in 1:num_folds) {
    f <- unique(data[, fold])[i]
    train <- data[which(data[, fold] == f), ]
    suppressMessages(
      train_fit <- logreg2ph(Y_unval = Y_unval, Y = Y, X_unval = X_unval, X = X, Z = Z,
                                        Bspline = Bspline, data = train,
                                        theta_pred = theta_pred, gamma_pred = gamma_pred,
                                        noSE = TRUE, TOL = TOL, MAX_ITER = MAX_ITER)
    )
    status[i] <- train_fit$converged
    msg[i] <- train_fit$converged_msg

    if (train_fit$converged) {
      train_theta <- train_fit$coeff$coeff
      train_gamma <- train_fit$outcome_err_coeff$coeff
      train_p <- train_fit$Bspline_coeff
      train_x <- data.frame(train[train[, Validated] == 1, X])
      train_x <- data.frame(train_x[order(train_x[, 1]), ])
      colnames(train_x) <- X
      train_x <- cbind(k = 1:nrow(train_x), train_x)
      train_p <- merge(train_x, train_p)

      test <- data[which(data[, fold] != f), ]
      test_x <- data.frame(test[test[, Validated] == 1, X])
      test_x <- data.frame(test_x[order(test_x[, 1]), ])
      colnames(test_x) <- X
      test_x <- cbind(k_ = 1:nrow(test_x), test_x)
      test_p <- matrix(data = NA, nrow = nrow(test_x), ncol = length(Bspline))

      for (i in 1:nrow(test_x)) {
        x_ <- test_x[i, X]
        bf <- suppressWarnings(expr = max(which(train_x[, X] <= x_)))
        af <- suppressWarnings(expr = min(which(train_x[, X] >= x_)))
        if (bf == -Inf) { bf <- af }
        if (af == Inf) { af <- bf }

        # x values immediately before/after
        x0 <- train_p[bf, X]
        x1 <- train_p[af, X]

        # B-spline coefficients immediately before/after
        p0 <- train_p[bf, -c(1:(1 + length(X)))]
        p1 <- train_p[af, -c(1:(1 + length(X)))]

        if (x1 == x0) {
          test_p[i, ] <- unlist(p0)
        } else {
          test_p[i, ] <- unlist((p0 * (x1 - x_) + p1 * (x_ - x0)) / (x1 - x0))
        }
      }

      # Recale columns of test_p to sum to 1
      denom <- colSums(test_p)
      denom[denom == 0] <- 1 # Avoid NaN error due to dividing by 0
      re_test_p <- t(t(test_p) / denom)

      # Construct complete dataset
      cd <- complete_data(Y_unval = "Ystar", Y = "Y", X_unval = "Xbstar", X = "Xb", Z = "Xa",
                          Validated = "V", Bspline = Bspline, data = test)
      # Calculate log-likelihood -------------------------------------------
      # ll <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
      #                            Y_unval = Y_unval, Y = Y, X_unval = X_unval, X = X, Z = Z,
      #                            Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
      #                            theta = train_theta, gamma = train_gamma, p = test_p)
      # fold_ll <- append(fold_ll, ll)

      ll_f <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
                                 Y_unval = Y_unval, Y = Y, X_unval = X_unval, X = X, Z = Z,
                                 Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
                                 theta = train_theta, gamma = train_gamma, p = re_test_p)
      ll[i] <- ll_f
    } else {

    }
  }
  return(list(loglik = ll, status = status, msg = msg))
}

