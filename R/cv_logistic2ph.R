#' Cross-validation log-likelihood prediction for \code{logistic2ph}
#'
#' Performs cross-validation to calculate the average predicted log likelihood for the \code{logistic2ph} function. This function can be used to select the B-spline basis that yields the largest average predicted log likelihood. See pacakge vigenette for code examples.
#'
#' @param y_unval Column name of the error-prone or unvalidated binary outcome. This argument is optional. If \code{y_unval = NULL} (the default), \code{y} is treated as error-free.
#' @param y Column name that stores the validated value of \code{y_unval} in the second phase. Subjects with missing values of \code{y} are considered as those not selected in the second phase. This argument is required.
#' @param x_unval Specifies the columns of the error-prone covariates. This argument is required.
#' @param x Specifies the columns that store the validated values of \code{x_unval} in the second phase. Subjects with missing values of \code{x} are considered as those not selected in the second phase. This argument is required.
#' @param z Specifies the columns of the accurately measured covariates. Subjects with missing values of \code{z} are omitted from the analysis. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param nfolds Specifies the number of cross-validation folds. The default value is \code{5}. Although \code{nfolds} can be as large as the sample size (leave-one-out cross-validation), it is not recommended for large datasets. The smallest value allowable is \code{3}.
#' @param max_iter Specifies the maximum number of iterations in the EM algorithm. The default number is \code{2000}. This argument is optional.
#' @param tol Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#'
#' @details
#' \code{cv_logistic2ph} gives log-likelihood prediction for models and data like those in \code{logistic2ph}. Therefore, the arguments of \code{cv_logistic2ph} is analogous to that of \code{logistic2ph}.
#'
#' @return
#' `cv_logistic2ph()` returns a list that includes the following components:
#' \item{avg_pred_loglike}{Stores the average predicted log likelihood.}
#' \item{pred_loglike}{Stores the predicted log likelihood in each fold.}
#' \item{converge}{Stores the convergence status of the EM algorithm in each run.}
#' @importFrom Rcpp evalCpp
#' @importFrom stats pchisq
#'
#' @export
cv_logistic2ph <- function(y_unval = NULL, y = NULL, x_unval = NULL, x = NULL, z = NULL,
                          data, nfolds = 5, tol = 1E-4, max_iter = 1000,
                           verbose = FALSE) {

  # variable name change
  Y_unval = y_unval; Y = y ; X_unval = x_unval; X = x; Z = z
  Bspline =  attr(data, "bs_name"); TOL = tol; MAX_ITER = max_iter

  if (nfolds >= 3) {
    if (verbose) {
      print(paste0(nfolds, "-folds cross-validation will be performed."))
    }
  } else {
    stop("nfolds needs to be greater than or equal to 3!")
  }

  # Create vector of predictors for each model ----------------------
  theta_pred <- c(X, Z)
  gamma_pred <- c(X_unval, Y, X, Z)
  # ----------------------  Create vector of predictors for each model

  # Define the fold column ------------------------------------------
  data$fold = sample(x = 1:nfolds,
                     size = nrow(data),
                     replace = TRUE)
  # ------------------------------------------ Define the fold column

  # Loop over the folds ---------------------------------------------
  status <- rep(TRUE, nfolds)
  msg <- rep("", nfolds)
  ll <- rep(NA, nfolds)
  for (j in 1:nfolds) {
    f <- unique(data[, "fold"])[j]
    train <- data[which(data[, "fold"] == f), ]
    suppressMessages(
      train_fit <- logistic2ph_all(y_unval = Y_unval, y = Y, x_unval = X_unval,
                                   x = X, z = Z, b_spline = Bspline, data = train,
                                   se = FALSE, tol = TOL, max_iter = MAX_ITER)
    )
    status[j] <- train_fit$converge

    if (train_fit$converge) {
      train_theta <- train_fit$coefficients$Estimate
      train_gamma <- train_fit$outcome_err_coefficients$Estimate
      train_p <- train_fit$Bspline_coeff ## k, x, and coefficients
      #train_x <- train_fit$Bspline_coefficients[, c(1:2)] ## k and x vals
      # train_x <- data.frame(train[which(!is.na(train[, X])), X])
      # train_x <- data.frame(train_x[order(train_x[, 1]), ])
      # colnames(train_x) <- X
      # train_x <- cbind(k = 1:nrow(train_x), train_x)
      # train_p <- merge(train_x, train_p)

      test <- data[which(data[, "fold"] != f), ]
      test_x <- unique(data.frame(test[which(!is.na(test[, X])), X]))
      test_x <- data.frame(test_x[order(test_x[, 1]), ])
      colnames(test_x) <- X
      test_x <- cbind(k_ = 1:nrow(test_x), test_x)
      test_p <- matrix(data = NA,
                       nrow = nrow(test_x),
                       ncol = length(Bspline))

      for (i in 1:nrow(test_x)) {
        x_ <- test_x[i, X]
        bf <- suppressWarnings(expr = max(which(train_p[, X] <= x_)))
        af <- suppressWarnings(expr = min(which(train_p[, X] >= x_)))
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
      test$V <- as.numeric(!is.na(test[, X])) ## define validation indicator
      cd <- complete_data(Y_unval = Y_unval,
                          Y = Y,
                          X_unval = X_unval,
                          X = X,
                          Z = Z,
                          Validated = "V",
                          Bspline = Bspline,
                          data = test)

      # Calculate log-likelihood -------------------------------------------
      ll_f <- observed_data_loglik(N = nrow(test),
                                   n = sum(!is.na(test[, X])),
                                   Y_unval = Y_unval,
                                   Y = Y,
                                   X_unval = X_unval,
                                   X = X,
                                   Z = Z,
                                   Bspline = Bspline,
                                   comp_dat_all = cd,
                                   theta_pred = theta_pred,
                                   gamma_pred = gamma_pred,
                                   theta = train_theta,
                                   gamma = train_gamma,
                                   p = re_test_p)
      ll[j] <- ll_f
    } else {
      ll[j] <- NA
    }
  }
  # --------------------------------------------- Loop over the folds

  # Return the results ----------------------------------------------
  avg_pred_loglik = mean(ll, na.rm = TRUE)
  res_final = list(avg_pred_loglik = avg_pred_loglik,
                   pred_loglik = ll,
                   converge = status)
  return(res_final)
  # ---------------------------------------------- Return the results
}
