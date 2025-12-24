utils::globalVariables(c("..X", ".I", "k", ":=")) # Silence CHECK issues with data.table

#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2021). See pacakge vigenette for code examples.
#'
#' @param y_unval Column name of the error-prone or unvalidated binary outcome. This argument is optional. If \code{y_unval = NULL} (the default), \code{y} is treated as error-free.
#' @param y Column name that stores the validated value of \code{y_unval} in the second phase. Subjects with missing values of \code{y} are considered as those not selected in the second phase. This argument is required.
#' @param x_unval Specifies the columns of the error-prone covariates. This argument is required.
#' @param x Specifies the columns that store the validated values of \code{x_unval} in the second phase. Subjects with missing values of \code{x} are considered as those not selected in the second phase. This argument is required.
#' @param z Specifies the columns of the accurately measured covariates. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param hn_scale Specifies the scale of the perturbation constant in the variance estimation. For example, if \code{hn_scale = 0.5}, then the perturbation constant is \eqn{0.5n^{-1/2}}, where \eqn{n} is the first-phase sample size. The default value is \code{1}. This argument is optional.
#' @param max_iter Maximum number of iterations in the EM algorithm. The default number is \code{1000}. This argument is optional.
#' @param tol Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param se If \code{FALSE}, then the variances of the parameter estimators will not be estimated. The default value is \code{TRUE}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#'
#' @details
#' Models for \code{logistic2ph()} are specified through the arguments. The dataset input should at least contain columns for unvalidated error-prone outcome, validated error-prone outcome,
#' unvalidated error-prone covariate(s), validated error-prone covariate(s), and B-spline basis. B-spline basis can be generated from \code{splines::bs()} function, with argument \code{x}
#' being the unvalidated error-prone covariate(s). See vignette for options in tuning the B-spline basis.
#'
#' @return
#' `logistic2ph()` returns an object of class `"logistic2ph"`. The function `coef()` is used to obtain the coefficients of the fitted model. The function `summary()` is used to obtain and print a summary of results.
#'
#' An object of class `"logistic2ph"` is a list containing at least the following components:
#' \item{call}{the matched call.}
#' \item{coefficients}{A named vector of the logistic regression coefficient estimates.}
#' \item{covariance}{The covariance matrix of the logistic regression coefficient estimates.}
#' \item{converge}{In parameter estimation, if the EM algorithm converges, then \code{converge = TRUE}. Otherwise, \code{converge = FALSE}.}
#' \item{converge_cov}{In variance estimation, if the EM algorithm converges, then \code{converge_cov = TRUE}. Otherwise, \code{converge_cov = FALSE}.}
#'
#' @examples
#' \dontrun{
#' # Regression model: ADE ~ CD4 + Prior_ART. ADE and CD4 are partially validated.
#' data("mock.vccc")
#' sn <- 20
#' data.logistic <- spline2ph(x = "CD4_unval", size = 20, degree = 3,
#'                            data = mock.vccc, group = "Prior_ART",
#'                            split_group = TRUE)
#' res_logistic <- logistic2ph(y = "ADE_val", y_unval = "ADE_unval",
#'                            x = "CD4_val", x_unval = "CD4_unval",
#'                             z = "Prior_ART", data = data.logistic,
#'                             hn_scale = 1/2, se = TRUE, tol = 1e-04,
#'                             max_iter = 1000, verbose = FALSE)
#' }
#'
#' @references
#' Lotspeich, S. C., Shepherd, B. E., Amorim, G. G. C., Shaw, P. A., & Tao, R. (2021). Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort. *Biometrics, biom.13512.* https://doi.org/10.1111/biom.13512
#'
#' @importFrom stats as.formula
#' @importFrom stats glm
#' @importFrom data.table data.table
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table setorderv
#' @importFrom data.table rbindlist
#' @importFrom data.table CJ
#' @importFrom data.table copy
#' @importFrom data.table ":="
#'
#' @export

logistic2ph <- function(
    y_unval  = NULL,
    y        = NULL,
    x_unval  = NULL,
    x        = NULL,
    z        = NULL,
    data     = NULL,
    hn_scale = 1,
    se       = TRUE,
    tol      = 1E-4,
    max_iter = 1000,
    verbose  = FALSE)
{
  # Store the function call
  model_call <- match.call()

  # variable name change
  Y_unval = y_unval; Y = y ; X_unval = x_unval; X = x; Z = z
  Bspline = attr(data, "bs_name"); noSE = !se; TOL = tol; MAX_ITER = max_iter

  if (missing(data))    stop("No dataset is provided!")
  if (missing(Bspline)) stop("The B-spline basis is not specified!")
  if (missing(Y))       stop("The accurately measured response Y is not specified!")
  if (xor(missing(X), missing(X_unval)))
    stop("If X_unval and X are NULL, all predictors are assumed to be error-free. You must define both variables or neither!")

  # Convert data to data.frame to avoid conflicts with tbl
  data <- data.frame(data)

  # Save sample size
  N <- nrow(data)

  # Calculate the validated subjects
  Validated <- complete.cases(data[,c(Y,X)])

  # n is how many validated subjects there are in data
  n <- sum(Validated)

  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(Validated), decreasing = TRUE), ]

  # Create vector of predictors for each model ----------------------
  theta_pred <- c(X, Z)
  gamma_pred <- c(X_unval, Y, X, Z)
  pred       <- unique(c(theta_pred, gamma_pred))

  # Determine error setting -----------------------------------------
  ## If unvalidated outcome was left blank, assume error-free -------
  errorsY    <- !is.null(Y_unval)
  if (!errorsY && any(is.na(data[, Y])))
    stop("If Y_unval is NULL, the outcome is assumed to be error-free, but Y has missing values.")

  # Add the B spline basis ------------------------------------------
  sn <- ncol(data[, Bspline])
  if(0 %in% colSums(data[c(1:n), Bspline]))
  {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

    res_coefficients <- data.frame(Estimate = rep(NA, length(theta_pred) + 1),
                                   SE = NA,
                                   Statistic = NA,
                                   pvalue = NA)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    return(list(coefficients = res_coefficients,
                covariance   = NA,
                converge     = FALSE,
                converge_cov = NA))
  }

  # Create "complete" data ------------------------------------------
  ## Save distinct X ------------------------------------------------
  setDT(data)

  # --- Step 1: Distinct and ordered x_obs ----------------------------
  if(verbose) message("Create distinct and ordered x_obs")
  x_obs <- unique(data[1:n, ..X])
  setkeyv(x_obs, X)            # sets key and orders in place
  m <- nrow(x_obs)

  # --- Step 2: Replicate x_obs (stacked)
  if(verbose) message("Replicate x_obs (stacked)")
  # Using .SD for in-place repetition avoids copies
  x_obs_stacked <- x_obs[rep(seq_len(m), times = (N - n))]
  # The key on X remains valid, no need to re-sort
  # (if order matters strictly, uncomment next line)
  setorderv(x_obs_stacked, X)

  # --- Step 3: comp_dat_val ------------------------------------------
  # Static component (first n observations)
  if(verbose) message("Locate static component (first n observations)")
  comp_dat_val <- data[1:n, c(Y_unval, X_unval, Z, Bspline, X, Y), with = FALSE]

  # Prepare keyed x_obs with index k
  x_obs[, k := .I]             # fast in-place index column

  # Set keys for fast join
  setkeyv(comp_dat_val, X)     # join key
  # Join directly by key — fastest possible join
  comp_dat_val <- x_obs[comp_dat_val, nomatch = 0L]

  # Select relevant columns
  comp_dat_val <- comp_dat_val[, c(Y_unval, pred, Bspline, "k"), with = FALSE]

  # Convert to matrix only at the end (copy unavoidable)
  comp_dat_val <- as.matrix(comp_dat_val)

  # --- Step 4: comp_dat_unval ----------------------------------------
  if(verbose) message("Efficient replication: (N - n) rows per x_obs")
  # Using CJ() cross join to generate index pairs without materializing full cross product in R memory
  idx <- CJ(k = seq_len(m), i = seq(from = n + 1, to = N))
  # Join these index pairs to data and x_obs
  comp_dat_unval <- data.table(
    data[rep((n + 1):N, times = m), c(Y_unval, setdiff(pred, X), Bspline), with = FALSE],
    x_obs_stacked,
    k = rep(seq_len(m), each = (N - n))
  )

  if (errorsY)
  {
    ### Create two versions of the complete data
    ### Create two independent copies of comp_dat_unval
    comp_dat_y0 <- copy(comp_dat_unval)
    comp_dat_y1 <- copy(comp_dat_unval)

    ### One with Y = 0
    comp_dat_y0[, (Y) := 0]

    ### One with Y = 1
    comp_dat_y1[, (Y) := 1]

    ### Combine vertically
    comp_dat_unval <- rbindlist(list(comp_dat_y0, comp_dat_y1))

    ### Convert to matrix if needed
    comp_dat_unval <- as.matrix(comp_dat_unval)
  }
  # Select columns (note: if already matrix, indexing by column names may fail)
  comp_dat_unval <- comp_dat_unval[, c(Y_unval, pred, Bspline, "k"), drop = FALSE]

  # Combine with comp_dat_val
  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)

  # Initialize B-spline coefficients {p_kj}  -----------------------
  ## Numerators sum B(Xi*) over k = 1,...,m ------------------------
  ## Save as p_val_num for updates (contributions don't change) ----
  p_val_num <- rowsum(x = comp_dat_val[, Bspline],
                      group = comp_dat_val[, "k"],
                      reorder = TRUE)
  prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))

  # Create model formulas
  theta_formula    <- as.formula(paste0(Y, "~", paste(theta_pred, collapse = "+")))
  theta_design_mat <- as.matrix(cbind(int = 1, comp_dat_all[, theta_pred]))
  gamma_formula    <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
  gamma_design_mat <- as.matrix(cbind(int = 1, comp_dat_all[, gamma_pred]))

  # Initialize parameter values -------------------------------------
  prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
  prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)

  CONVERGED <- FALSE
  it <- 1

  # Estimate theta using EM -------------------------------------------
  if(verbose) message("Beginning EM loop")

  while(it <= MAX_ITER && !CONVERGED)
  {
    if(verbose) message("Iteration: ", it, "  Time: ", Sys.time())

    # ===================================================================
    # E-Step (C++)
    if (!is.matrix(comp_dat_unval))
      comp_dat_unval <- as.matrix(comp_dat_unval)

    Bspline_mat <- comp_dat_unval[, Bspline, drop = FALSE]
    storage.mode(Bspline_mat) <- "double"

    Y_vec <- comp_dat_unval[, Y]
    Y_unval_vec <- if (errorsY) comp_dat_unval[, Y_unval] else rep(0, nrow(comp_dat_unval))
    storage.mode(Y_vec) <- "double"
    storage.mode(Y_unval_vec) <- "double"

    res_C <- logistic2ph_estep(
      theta_design_mat,
      gamma_design_mat,
      prev_theta,
      prev_gamma,
      prev_p,
      Bspline_mat,
      Y_vec,
      Y_unval_vec,
      n, N, m,
      errorsY
    )

    w_t <- res_C$w_t
    u_t <- res_C$u_t

    # ===================================================================
    # M-Step (C++)
    ## Lengthen w_t by prepending n ones (for validated subjects)
    w_t <- .lengthenWT(w_t, n)

    mstep_result <- logistic2ph_mstep(
      theta_design_mat,
      gamma_design_mat,
      comp_dat_all[, Y],
      comp_dat_all[, Y_unval],
      w_t,
      u_t,
      p_val_num,
      prev_theta,
      prev_gamma,
      prev_p,
      n, N, m,
      errorsY,
      TOL
    )

    # Extract results
    new_theta <- mstep_result$theta
    new_gamma <- mstep_result$gamma
    new_p     <- mstep_result$p

    # ===================================================================
    # Fallback to R glm() if C++ failed (singular Hessian)
    if (!mstep_result$theta_success) {
      if(verbose) warning("C++ M-step failed for theta, falling back to glm")
      suppressWarnings(
        new_theta <- matrix(
          glm(formula = theta_formula,
              family = "binomial",
              data = data.frame(comp_dat_all),
              weights = w_t)$coefficients,
          ncol = 1
        )
      )
    }

    if (errorsY && !mstep_result$gamma_success) {
      if(verbose) warning("C++ M-step failed for gamma, falling back to glm")
      suppressWarnings(
        new_gamma <- matrix(
          glm(formula = gamma_formula,
              family = "binomial",
              data = data.frame(comp_dat_all),
              weights = w_t)$coefficients,
          ncol = 1
        )
      )
    }

    # ===================================================================
    # Check convergence
    if (mstep_result$theta_conv && mstep_result$gamma_conv && mstep_result$p_conv)
      CONVERGED <- TRUE

    # ===================================================================
    # Error checking for singular systems
    theta_conv <- abs(new_theta - prev_theta) < TOL
    if(any(is.na(theta_conv)))
    {
      bad <- paste(all.vars(theta_formula)[is.na(theta_conv)], collapse=", ")
      stop(paste("System is singular. Variables possibly at fault: ", bad))
    }

    if (errorsY)
    {
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      if(any(is.na(gamma_conv)))
      {
        bad <- paste(all.vars(gamma_formula)[is.na(gamma_conv)], collapse=", ")
        stop(paste("System is singular. Variables possibly at fault: ", bad))
      }
    }

    # ===================================================================
    # Update for next iteration
    it <- it + 1
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p
  }

  # =====================================================================
  # Post-convergence processing
  # =====================================================================
  rownames(new_theta) <- c("Intercept", theta_pred)
  if (errorsY) rownames(new_gamma) <- c("Intercept", gamma_pred)

  if(!CONVERGED) {
    if(verbose && it > MAX_ITER) message("MAX_ITER reached")

    res_coefficients <- data.frame(
      Estimate  = rep(NA, length(new_theta)),
      SE        = NA,
      Statistic = NA,
      pvalue    = NA
    )
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    return(list(
      coefficients = res_coefficients,
      covariance   = NA,
      converge     = FALSE,
      converge_cov = NA
    ))
  }

  if(verbose) message("EM algorithm converged in ", it - 1, " iterations")

  # If computation of SE is not requested, return results
  if(noSE)
  {
    res_coefficients <- data.frame(Estimate  = new_theta,
                                   SE        = NA,
                                   Statistic = NA,
                                   pvalue    = NA)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    return(list(coefficients = res_coefficients,
                covariance = NA,
                converge = CONVERGED,
                converge_cov = NA))
  }

  # Estimate Cov(theta) using profile likelihood -------------------------
  h_N <- hn_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

  ## Calculate pl(theta) -------------------------------------------------
  od_loglik_theta <- observed_data_loglik(
    N = N,
    n = n,
    Y_unval = Y_unval,
    Y = Y,
    X_unval = X_unval,
    X = X,
    Z = Z,
    Bspline = Bspline,
    comp_dat_all = comp_dat_all,
    theta_pred = theta_pred,
    gamma_pred = gamma_pred,
    theta = new_theta,
    gamma = new_gamma,
    p = new_p)

  I_theta <- matrix(data = od_loglik_theta,
                    nrow = nrow(new_theta),
                    ncol = nrow(new_theta))

  single_pert_theta <- sapply(X = seq(1, ncol(I_theta)),
    FUN = pl_theta,
    theta = new_theta,
    h_N = h_N,
    n = n,
    N = N,
    Y_unval = Y_unval,
    Y = Y,
    X_unval = X_unval,
    X,
    Z = Z,
    Bspline = Bspline,
    comp_dat_all = comp_dat_all,
    theta_pred = theta_pred,
    gamma_pred = gamma_pred,
    gamma0 = new_gamma,
    p0 = new_p,
    p_val_num = p_val_num,
    TOL = TOL,
    MAX_ITER = MAX_ITER)

  if (any(is.na(single_pert_theta)))
  {
    I_theta <- matrix(data = NA,
                      nrow = nrow(new_theta),
                      ncol = nrow(new_theta))
    SE_CONVERGED <- FALSE
  } else
  {
    spt_wide <- matrix(data = rep(c(single_pert_theta),
                           times = ncol(I_theta)),
                       ncol = ncol(I_theta),
                       byrow = FALSE)
    #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
    I_theta <- I_theta - spt_wide - t(spt_wide)
    SE_CONVERGED <- TRUE
  }

  for (c in 1:ncol(I_theta))
  {
    pert_theta <- new_theta
    pert_theta[c] <- pert_theta[c] + h_N
    double_pert_theta <- sapply(X = seq(c, ncol(I_theta)),
      FUN = pl_theta,
      theta = pert_theta,
      h_N = h_N,
      n = n,
      N = N,
      Y_unval = Y_unval,
      Y = Y,
      X_unval = X_unval,
      X,
      Z = Z,
      Bspline = Bspline,
      comp_dat_all = comp_dat_all,
      theta_pred = theta_pred,
      gamma_pred = gamma_pred,
      gamma0 = new_gamma,
      p0 = new_p,
      p_val_num = p_val_num,
      MAX_ITER = MAX_ITER,
      TOL = TOL)
    dpt <- matrix(data = 0,
                  nrow = nrow(I_theta),
                  ncol = ncol(I_theta))
    dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
    if(c < ncol(I_theta))
    {
      ## And fill the others in on the cth row/ column
      dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
    }

    I_theta <- I_theta + dpt
  }
  I_theta <- h_N ^ (- 2) * I_theta
  cov_theta <- tryCatch(expr = - solve(I_theta),
    error = function(err) { matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta)) }
    )
  # ------------------------- Estimate Cov(theta) using profile likelihood

  se_theta <- tryCatch(
    expr    = sqrt(diag(cov_theta)),
    warning = function(w) { matrix(NA, nrow = nrow(prev_theta)) }
  )
  SE_CONVERGED <- !any(is.na(se_theta))

  res_coefficients           <- data.frame(Estimate = new_theta, SE = se_theta)
  res_coefficients$Statistic <- res_coefficients$Estimate / res_coefficients$SE
  res_coefficients$pvalue    <- 1 - pchisq(res_coefficients$Statistic ^ 2, df = 1)
  colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

  coef_return                <- as.numeric(new_theta)
  names(coef_return)         <- rownames(res_coefficients)

  logistic2ph_class(list(
    call         = model_call,
    coefficients = coef_return,
    covariance   = cov_theta,
    converge     = CONVERGED,
    converge_cov = SE_CONVERGED))
}
