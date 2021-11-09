#' Cross-validated observed-data log-likelihood
#' This function returns the value of the observed-data log-likelihood based on cross-validation.
#'
#' @param fold Column name with the assigned fold for cross-validation.
#' @param mod_Y_unvalidated Column name with the unvalidated outcome. If \code{mod_Y_unvalidated} is null, the outcome is assumed to be error-free.
#' @param mod_Y_validated Column name with the validated outcome.
#' @param mod_X_unvalidated Column name(s) with the unvalidated predictors.  If \code{mod_X_unvalidated} and \code{mod_X_validated} are \code{null}, all precictors are assumed to be error-free.
#' @param mod_X_validated Column name(s) with the validated predictors. If \code{mod_X_unvalidated} and \code{mod_X_validated} are \code{null}, all precictors are assumed to be error-free.
#' @param true_covariates (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{mod_Y_unvalidated}, \code{mod_Y_validated}, \code{mod_X_unvalidated}, \code{mod_X_validated}, \code{true_covariates}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return scalar value of the function
#' @export

cv_observed_data_loglik <- function(fold, mod_Y_unvalidated = NULL, mod_Y_validated = NULL, mod_X_unvalidated = NULL, mod_X_validated = NULL, true_covariates = NULL,
                        Validated = NULL, Bspline = NULL, data, theta_pred = NULL, gamma_pred = NULL,
                        TOL = 1E-4, MAX_ITER = 1000) {
  if (is.null(theta_pred)) { theta_pred <- c(mod_X_validated, true_covariates) }
  if (is.null(gamma_pred) & !is.null(mod_Y_unvalidated)) { gamma_pred <- c(mod_X_unvalidated, mod_Y_validated, mod_X_validated, true_covariates) }

  num_folds <- length(unique(data[, fold]))
  status <- rep(TRUE, num_folds)
  msg <- rep("", num_folds)
  ll <- rep(NA, num_folds)
  #fold_ll <- re_fold_ll <- vector()
  for (i in 1:num_folds) {
    f <- unique(data[, fold])[i]
    train <- data[which(data[, fold] == f), ]
    suppressMessages(
      train_fit <- logreg2ph::logreg2ph(mod_Y_unvalidated = mod_Y_unvalidated, mod_Y_validated = mod_Y_validated, mod_X_unvalidated = mod_X_unvalidated, mod_X_validated = mod_X_validated, true_covariates = true_covariates,
                                        Validated = Validated, Bspline = Bspline, data = train,
                                        theta_pred = theta_pred, gamma_pred = gamma_pred,
                                        noSE = TRUE, TOL = TOL, MAX_ITER = MAX_ITER)
    )
    status[i] <- train_fit$converged
    msg[i] <- train_fit$converged_msg

    if (train_fit$converged) {
      train_theta <- train_fit$coeff$coeff
      train_gamma <- train_fit$outcome_err_coeff$coeff
      train_p <- train_fit$Bspline_coeff
      train_x <- data.frame(train[train[, Validated] == 1, mod_X_validated])
      train_x <- data.frame(train_x[order(train_x[, 1]), ])
      colnames(train_x) <- mod_X_validated
      train_x <- cbind(k = 1:nrow(train_x), train_x)
      train_p <- merge(train_x, train_p)

      test <- data[which(data[, fold] != f), ]
      test_x <- data.frame(test[test[, Validated] == 1, mod_X_validated])
      test_x <- data.frame(test_x[order(test_x[, 1]), ])
      colnames(test_x) <- mod_X_validated
      test_x <- cbind(k_ = 1:nrow(test_x), test_x)
      test_p <- matrix(data = NA, nrow = nrow(test_x), ncol = length(Bspline))

      for (i in 1:nrow(test_x)) {
        x_ <- test_x[i, mod_X_validated]
        bf <- suppressWarnings(expr = max(which(train_x[, mod_X_validated] <= x_)))
        af <- suppressWarnings(expr = min(which(train_x[, mod_X_validated] >= x_)))
        if (bf == -Inf) { bf <- af }
        if (af == Inf) { af <- bf }

        # x values immediately before/after
        x0 <- train_p[bf, mod_X_validated]
        x1 <- train_p[af, mod_X_validated]

        # B-spline coefficients immediately before/after
        p0 <- train_p[bf, -c(1:(1 + length(mod_X_validated)))]
        p1 <- train_p[af, -c(1:(1 + length(mod_X_validated)))]

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
      cd <- complete_data(mod_Y_unvalidated = "Ystar", mod_Y_validated = "Y", mod_X_unvalidated = "Xbstar", mod_X_validated = "Xb", true_covariates = "Xa",
                          Validated = "V", Bspline = colnames(B), data = test)
      # Calculate log-likelihood -------------------------------------------
      # ll <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
      #                            mod_Y_unvalidated = mod_Y_unvalidated, mod_Y_validated = mod_Y_validated, mod_X_unvalidated = mod_X_unvalidated, mod_X_validated = mod_X_validated, true_covariates = true_covariates,
      #                            Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
      #                            theta = train_theta, gamma = train_gamma, p = test_p)
      # fold_ll <- append(fold_ll, ll)

      ll_f <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
                                 mod_Y_unvalidated = mod_Y_unvalidated, mod_Y_validated = mod_Y_validated, mod_X_unvalidated = mod_X_unvalidated, mod_X_validated = mod_X_validated, true_covariates = true_covariates,
                                 Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
                                 theta = train_theta, gamma = train_gamma, p = re_test_p)
      ll[i] <- ll_f
    } else {

    }
  }
  return(list(loglik = ll, status = status, msg = msg))
}

