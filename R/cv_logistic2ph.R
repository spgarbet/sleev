#' Performs cross-validation to calculate the average predicted log likelihood for the \code{logistic2ph} function. This function can be used to select the B-spline basis that yields the largest average predicted log likelihood.
#'
#' @param Y_unval Column name of the error-prone or unvalidated binary outcome. This argument is optional. If \code{Y_unval = NULL} (the default), \code{Y} is treated as error-free.
#' @param Y Column name that stores the validated value of \code{Y_unval} in the second phase. Subjects with missing values of \code{Y} are considered as those not selected in the second phase. This argument is required.
#' @param X_unval Specifies the columns of the error-prone covariates. This argument is required.
#' @param X Specifies the columns that store the validated values of \code{X_unval} in the second phase. Subjects with missing values of \code{X} are considered as those not selected in the second phase. This argument is required.
#' @param Bspline Specifies the columns of the B-spline basis. Subjects with missing values of \code{Bspline} are omitted from the analysis. This argument is required.
#' @param Z Specifies the columns of the accurately measured covariates. Subjects with missing values of \code{Z} are omitted from the analysis. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param nfolds Specifies the number of cross-validation folds. The default value is \code{5}. Although \code{nfolds} can be as large as the sample size (leave-one-out cross-validation), it is not recommended for large datasets. The smallest value allowable is \code{3}.
#' @param MAX_ITER Specifies the maximum number of iterations in the EM algorithm. The default number is \code{2000}. This argument is optional.
#' @param TOL Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#' @return
#' \item{avg_pred_loglike}{Stores the average predicted log likelihood.}
#' \item{pred_loglike}{Stores the predicted log likelihood in each fold.}
#' \item{converge}{Stores the convergence status of the EM algorithm in each run.}
#' @importFrom Rcpp evalCpp
#' @importFrom stats pchisq
#'
#' @export
cv_logistic2ph <- function(Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL,
                           Bspline = NULL, data, nfolds = 5, TOL = 1E-4, MAX_ITER = 1000,
                           verbose = FALSE) {
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
  for (i in 1:nfolds) {
    f <- unique(data[, "fold"])[i]
    train <- data[which(data[, "fold"] == f), ]
    suppressMessages(
      train_fit <- logistic2ph_all(Y_unval = Y_unval, Y = Y, X_unval = X_unval,
                                   X = X, Z = Z, Bspline = Bspline, data = train,
                                   noSE = TRUE, TOL = TOL, MAX_ITER = MAX_ITER)
    )
    status[i] <- train_fit$converge

    if (train_fit$converge) {
      train_theta <- train_fit$coeff$coeff
      train_gamma <- train_fit$outcome_err_coeff$coeff
      train_p <- train_fit$Bspline_coeff
      train_x <- data.frame(train[which(!is.na(train[, X])), X])
      # train_x <- data.frame(train[train[, Validated] == 1, X])
      train_x <- data.frame(train_x[order(train_x[, 1]), ])
      colnames(train_x) <- X
      train_x <- cbind(k = 1:nrow(train_x), train_x)
      train_p <- merge(train_x, train_p)

      test <- data[which(data[, "fold"] != f), ]
      test_x <- data.frame(test[which(!is.na(test[, X])), X])
      # test_x <- data.frame(test[test[, Validated] == 1, X])
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
      ll[i] <- ll_f
    } else {
      ll[i] <- NA
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
