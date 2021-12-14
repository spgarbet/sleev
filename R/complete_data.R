#' @title
#' Create complete dataset
#'
#' @description
#' This function returns the complete dataset trying all validated values for each unvalidated subject.
#'
#' @param Y_unval Column names with the unvalidated outcome. If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y Column names with the validated outcome.
#' @param X_unval Column name(s) with the unvalidated predictors.  If \code{X_unval} and \code{X} are \code{null}, all predictors are assumed to be error-free.
#' @param X Column name(s) with the validated predictors. If \code{X_unval} and \code{X} are \code{null}, all predictors are assumed to be error-free.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y}, \code{X_unval}, \code{X}, \code{Z}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @return dataframe

complete_data <- function(Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL,
                          Validated = NULL, Bspline = NULL, data, theta_pred = NULL, gamma_pred = NULL) {
  N <- nrow(data)
  n <- sum(data[, Validated])
  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(data[, Validated]), decreasing = TRUE), ]
  data$id <- 1:N
  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsY <- errorsX <- TRUE
  if (is.null(Y_unval)) {errorsY <- FALSE}
  if (is.null(X_unval) & is.null(X)) {errorsX <- FALSE}

  if (is.null(theta_pred)) {
    theta_pred <- Z(X, Z)
  }

  if (is.null(gamma_pred) & errorsY) {
    gamma_pred <- Z(X_unval, Y, X, Z)
  }

  pred <- unique(Z(theta_pred, gamma_pred))

  if (errorsX & errorsY) {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, Z(X)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- Z(X)

    # Save static (X*,Y*,X,Y,Z) since they don't change ---------------
    comp_dat_val <- data[Z(1:n), Z("id", Y_unval, X_unval, Z, Bspline, X, Y)]
    comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
    comp_dat_val <- comp_dat_val[, Z("id", Y_unval, pred, Bspline, "k")]
    comp_dat_val <- data.matrix(comp_dat_val)

    # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
    # one row per x) --------------------------------------------------
    suppressWarnings(comp_dat_unval <- cbind(data[-Z(1:n), Z("id", Y_unval, setdiff(x = pred, y = Z(Y, X)), Bspline)],
                                             x_obs_stacked))
    comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
    comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
    colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y
    comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1),
                                        k = rep(rep(seq(1, m), each = (N - n)), times = 2)))
    comp_dat_unval <- comp_dat_unval[, Z("id", Y_unval, pred, Bspline, "k")]

    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  } else if (errorsX) {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, Z(X)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- Z(X)

    # Save static (X*,X,Y,Z) since they don't change ---------------
    comp_dat_val <- data[Z(1:n), Z(Y, pred, Bspline)]
    comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
    comp_dat_val <- comp_dat_val[, Z(Y, pred, Bspline, "k")]
    comp_dat_val <- data.matrix(comp_dat_val)

    # (m x n)xd vectors of each (one column per person, one row per x) --
    suppressWarnings(
      comp_dat_unval <- data.matrix(
        cbind(data[-Z(1:n), Z(Y, setdiff(x = pred, y = Z(X)), Bspline)],
              x_obs_stacked,
              k = rep(seq(1, m), each = (N - n)))
      )
    )
    comp_dat_unval <- comp_dat_unval[, Z(Y, pred, Bspline, "k")]

    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  } else if (errorsY) {
    # Save static (Y*,X,Y,Z) since they don't change ------------------
    comp_dat_val <- data.matrix(data[Z(1:n), Z(Y_unval, pred)])

    # Create duplicate rows of each person (one each for y = 0/1) -----
    comp_dat_unval <- data[-Z(1:n), Z(Y_unval, setdiff(x = pred, y = Z(Y)))]
    comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
    comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
    colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <-
      colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y
    comp_dat_unval <- data.matrix(rbind(comp_dat_y0, comp_dat_y1))

    # Stack complete data: --------------------------------------------
    ## n rows for the n subjects in Phase II (1 each) -----------------
    ## 2 * (N - n) for the (N - n) subjects from Phase I (2 each) -----
    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  }
  #comp_dat_all <- cbind(id = Z(1:n, rep(seq(n+1, N), each = n)), comp_dat_all)
  return(comp_dat_all)
}
