#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2021).
#'
#' @param Y_unval Column name of the error-prone or unvalidated continuous outcome. Subjects with missing values of \code{Y_unval} are omitted from the analysis. If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y Column name that stores the validated value of \code{Y_unval} in the second phase. Subjects with missing values of \code{Y} are considered as those not selected in the second phase. This argument is required.
#' @param X_unval Column name(s) with the unvalidated predictors.  If \code{X_unval} and \code{X} are \code{null}, all predictors are assumed to be error-free.
#' @param X Column name(s) with the validated predictors. If \code{X_unval} and \code{X} are `NULL`, all predictors are assumed to be error-free.
#' @param Z (Optional) Column name(s) with additional error-free covariates.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y}, \code{X_unval}, \code{X}, \code{Z}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) \code{"Zeros"} (non-informative starting values) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param hn_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is `hn_scale = 1`.
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations in the EM algorithm. The default number is \code{1000}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#' 
#' @return
#' \item{coeff}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{outcome_err_coeff}{dataframe with final coefficient estimates for the outcome error model.}
#' \item{Bspline_coeff}{dataframe with final B-spline coefficient estimates (where applicable).}
#' \item{vcov}{variance-covarianced matrix for \code{coeff} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' \item{od_loglik_at_conv}{value of the observed-data log-likelihood at convergence.}
#' 
#' @references
#' Lotspeich, S. C., Shepherd, B. E., Amorim, G. G. C., Shaw, P. A., & Tao, R. (2021). Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort. *Biometrics, biom.13512.* https://doi.org/10.1111/biom.13512
#' 
#' @importFrom stats as.formula
#' @importFrom stats glm
#' 
#' @examples
#'  set.seed(918)
#'  
#'  # Set sample sizes ----------------------------------------
#'  N <- 1000 # Phase-I = N
#'  n <- 250 # Phase-II/audit size = n
#'  
#'  # Generate true values Y, Xb, Xa --------------------------
#'  Xa <- rbinom(n = N, size = 1, prob = 0.25)
#'  Xb <- rbinom(n = N, size = 1, prob = 0.5)
#'  Y <- rbinom(n = N, size = 1,prob = (1 + exp(-(- 0.65 - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))
#'  
#'  # Generate error-prone Xb* from error model P(Xb*|Xb,Xa) --
#'  sensX <- specX <- 0.75
#'  delta0 <- - log(specX / (1 - specX))
#'  delta1 <- - delta0 - log((1 - sensX) / sensX)
#'  Xbstar <- rbinom(n = N, size = 1,
#'                   prob = (1 + exp(- (delta0 + delta1 * Xb + 0.5 * Xa))) ^ (- 1))
#'  
#'  # Generate error-prone Y* from error model P(Y*|Xb*,Y,Xb,Xa)
#'  sensY <- 0.95
#'  specY <- 0.90
#'  theta0 <- - log(specY / (1 - specY))
#'  theta1 <- - theta0 - log((1 - sensY) / sensY)
#'  Ystar <- rbinom(n = N, size = 1,
#'    prob = (1 + exp(- (theta0 - 0.2 * Xbstar + theta1 * Y - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))
#'  
#'  
#'  ## V is a TRUE/FALSE vector where TRUE = validated --------
#'  V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
#'  
#'  
#'  # Build dataset --------------------------------------------
#'  sdat <- cbind(Y, Xb, Ystar, Xbstar, Xa)
#'  # Make Phase-II variables Y, Xb NA for unaudited subjects ---
#'  sdat[!V, c("Y", "Xb")] <- NA
#'  
#'  # Fit models -----------------------------------------------
#'  ## Naive model -----------------------------------------
#'  naive <- glm(Ystar ~ Xbstar + Xa, family = "binomial", data = data.frame(sdat))
#'  
#'  
#'  ## Generalized raking ----------------------------------
#'  ### Influence function for logistic regression
#'  ### Taken from: https://github.com/T0ngChen/multiwave/blob/master/sim.r
#'  inf.fun <- function(fit) {
#'    dm <- model.matrix(fit)
#'    Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
#'    ## influence function
#'    infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
#'    infl
#'  }
#'  naive_infl <- inf.fun(naive) # error-prone influence functions based on naive model
#'  colnames(naive_infl) <- paste0("if", 1:3)
#'  
#'  # Add naive influence functions to sdat -----------------------------------------------
#'  sdat <- cbind(id = 1:N, sdat, naive_infl)
#'  
#'  ### Construct B-spline basis -------------------------------
#'  ### Since Xb* and Xa are both binary, reduces to indicators --
#'  nsieve <- 4
#'  B <- matrix(0, nrow = N, ncol = nsieve)
#'  B[which(Xa == 0 & Xbstar == 0), 1] <- 1
#'  B[which(Xa == 0 & Xbstar == 1), 2] <- 1
#'  B[which(Xa == 1 & Xbstar == 0), 3] <- 1
#'  B[which(Xa == 1 & Xbstar == 1), 4] <- 1
#'  colnames(B) <- paste0("bs", seq(1, nsieve))
#'  sdat <- cbind(sdat, B)
#'  smle <- logistic2ph(Y_unval = "Ystar",
#'    Y = "Y",
#'    X_unval = "Xbstar",
#'    X = "Xb",
#'    Z = "Xa",
#'    Bspline = colnames(B),
#'    data = sdat,
#'    noSE = FALSE,
#'    MAX_ITER = 1000,
#'    TOL = 1E-4)
#' @export

logistic2ph <- function(Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL,Bspline = NULL, data = NULL, theta_pred = NULL, gamma_pred = NULL,initial_lr_params = "Zeros", hn_scale = 1, noSE = FALSE, TOL = 1E-4, MAX_ITER = 1000, verbose = FALSE)
{
  if (missing(data)) 
  {
    stop("No dataset is provided!")
  }

  if (missing(Bspline)) 
  {
    stop("The B-spline basis is not specified!")
  } 

  if (missing(Y)) 
  {
    stop("The accurately measured response Y is not specified!")
  }
  
  if (missing(X)) 
  {
    stop("The validated covariates in the second-phase are not specified!")
  }
  
  if (length(data[,X_unval]) != length(data[,X])) 
  {
    stop("The number of columns in X_unval and X is different!")
  }

  N <- nrow(data)

  # Calculate the validated subjects
  Validated <- logical(N) # initialize a logical vector of length N
  for (i in 1:N)
  {
    Validated[i] <- !(is.na(data[i,X]) || is.na(data[i,Y]))
  }
  
  # n is how many validated subjects there are in data
  n <- sum(Validated)
  
  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(Validated), decreasing = TRUE), ]

  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsY <- errorsX <- TRUE
  if (is.null(Y_unval))
  {errorsY <- FALSE}

  if (is.null(X_unval) & is.null(X))
  {errorsX <- FALSE}

  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  # Add the B spline basis ------------------------------------------
  if (errorsX)
  {
    sn <- ncol(data[, Bspline])
    if(0 %in% colSums(data[c(1:n), Bspline]))
    {
      warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

      return(list(coeff = data.frame(coeff = NA, se = NA),
        outcome_err_coeff = data.frame(coeff = NA, se = NA),
        Bspline_coeff = NA,
        vcov = NA,
        converged = NA,
        se_converged = NA,
        converged_msg = "B-spline error",
        iterations = 0,
        od_loglik_at_conv = NA))
    }

  }

  # ------------------------------------------ Add the B spline basis

  if (is.null(theta_pred))
  {
    theta_pred <- c(X, Z)
    if (verbose)
    {
      message("Analysis model assumed main effects only.")
    }
  }


  if (is.null(gamma_pred) & errorsY)
  {
    gamma_pred <- c(X_unval, Y, X, Z)
    if (verbose)
    {
      message("Outcome error model assumed main effects only.")
    }
  }


  pred <- unique(c(theta_pred, gamma_pred))

  if (errorsX & errorsY)
  {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, c(X)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- c(X)

    # Save static (X*,Y*,X,Y,Z) since they don't change ---------------
    comp_dat_val <- data[c(1:n), c(Y_unval, X_unval, Z, Bspline, X, Y)]
    comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
    comp_dat_val <- comp_dat_val[, c(Y_unval, pred, Bspline, "k")]
    comp_dat_val <- data.matrix(comp_dat_val)

    # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
    # one row per x) --------------------------------------------------
    suppressWarnings(comp_dat_unval <- cbind(data[-c(1:n), c(Y_unval, setdiff(x = pred, y = c(Y, X)), Bspline)],
     x_obs_stacked))
    comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
    comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
    colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y
    comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1),
      k = rep(rep(seq(1, m), each = (N - n)), times = 2)))
    comp_dat_unval <- comp_dat_unval[, c(Y_unval, pred, Bspline, "k")]

    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
    prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))
  }
  else if (errorsX)
  {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, c(X)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- c(X)

    # Save static (X*,X,Y,Z) since they don't change ---------------
    comp_dat_val <- data[c(1:n), c(Y, pred, Bspline)]
    comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
    comp_dat_val <- comp_dat_val[, c(Y, pred, Bspline, "k")]
    comp_dat_val <- data.matrix(comp_dat_val)

    # (m x n)xd vectors of each (one column per person, one row per x) --
    suppressWarnings(
      comp_dat_unval <- data.matrix(
        cbind(data[-c(1:n), c(Y, setdiff(x = pred, y = c(X)), Bspline)],
          x_obs_stacked,
          k = rep(seq(1, m), each = (N - n)))
        )
      )
    comp_dat_unval <- comp_dat_unval[, c(Y, pred, Bspline, "k")]

    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
    prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))
  }
  else if (errorsY)
  {
    # Save static (Y*,X,Y,Z) since they don't change ------------------
    comp_dat_val <- data.matrix(data[c(1:n), c(Y_unval, pred)])

    # Create duplicate rows of each person (one each for y = 0/1) -----
    comp_dat_unval <- data[-c(1:n), c(Y_unval, setdiff(x = pred, y = c(Y)))]
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


  theta_formula <- as.formula(paste0(Y, "~", paste(theta_pred, collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[, theta_pred])

  if (errorsY)
  {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- cbind(int = 1, comp_dat_all[, gamma_pred])
  }


  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zeros", "Complete-data")))
  {
    if (verbose)
    {
      message("Invalid starting values provided. Non-informative zeros assumed.")
    }
    initial_lr_params <- "Zeros"
  }


  if(initial_lr_params == "Zeros")
  {
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    if (errorsY)
    {
      prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
    }

    } else if(initial_lr_params == "Complete-data")
    {
      prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficients, ncol = 1)
      if (errorsY)
      {
        prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficient, ncol = 1)
      }

    }


    CONVERGED <- FALSE
    CONVERGED_MSG <- "Unknown"
    it <- 1

  # pre-allocate memory for loop variables
  mus_theta <- vector("numeric", nrow(theta_design_mat) * ncol(prev_theta))
  mus_gamma <- vector("numeric", nrow(gamma_design_mat) * ncol(prev_gamma))

  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric((theta_design_mat[-c(1:n), ] %*% prev_theta))
    pY_X <- 1 / (1 + exp(- mu_theta))
    I_y0 <- comp_dat_unval[, Y] == 0
    pY_X[I_y0] <- 1 - pY_X[I_y0]

    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (errorsY)
    {
      pYstar <- 1 / (1 + exp(- as.numeric((gamma_design_mat[-c(1:n), ] %*% prev_gamma))))
      I_ystar0 <- comp_dat_unval[, Y_unval] == 0
      pYstar[I_ystar0] <- 1 - pYstar[I_ystar0]

    } #else {
      #pYstar <- matrix(1, nrow = nrow(comp_dat_unval))
    #}

    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    # cpppX <- pXCalc(comp_dat_unval, errorsX, errorsY, prev_p, rep(seq(1, m), each = (N - n))-1, match(Bspline, colnames(comp_dat_unval))-1, seq(1, nrow(comp_dat_unval))-1)
    if (errorsX & errorsY)
    {
      ### p_kj ------------------------------------------------------
      ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
      ### multiply by the B-spline terms
      pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_unval[, Bspline]
      ### ---------------------------------------------------------- p_kj
    }
    else if (errorsX)
    {
      ### p_kj ----------------------------------------------------------
      ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
      ### multiply by the B-spline terms
      pX <- prev_p[rep(seq(1, m), each = (N - n)), ]  * comp_dat_unval[, Bspline]
      ### ---------------------------------------------------------- p_kj
    } #else if (errorsY) {
      #pX <- rep(1, times = nrow(comp_dat_unval))
    #}

    # tic("cond exp cpp")
    # condExp <- conditionalExpectations(errorsX, errorsY, pX, pY_X, pYstar, N-n, m)
    # CPPw_t <- condExp[["w_t"]]
    # CPPu_t <- condExp[["u_t"]]
    # CPPpsi_t <- condExp[["psi_t"]]
    # toc()

    # R is faster than cpp

    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    if (errorsY & errorsX)
    {
      ### P(Y|X,Z)P(Y*|X*,Y,X,Z)p_kjB(X*) -----------------------------
      psi_num <- c(pY_X * pYstar) * pX
      ### Update denominator ------------------------------------------
      #### Sum up all rows per id (e.g. sum over xk/y) ----------------
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2 * m))
      #### Then sum over the sn splines -------------------------------
      psi_denom <- rowSums(psi_denom)
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      ### by summing across the splines/ columns of psi_t -------------
      w_t <- rowSums(psi_t)
      ### Update the u_kji for unvalidated subjects ------------------
      ### by summing over Y = 0/1 w/i each i, k ----------------------
      ### add top half of psi_t (y = 0) to bottom half (y = 1) -------
      u_t <- psi_t[c(1:(m * (N - n))), ] + psi_t[- c(1:(m * (N - n))), ]
      } else if (errorsX)
      {
      ### P(Y|X,Z)p_kjB(X*) -------------------------------------------
      psi_num <- c(pY_X) * pX
      ### Update denominator ------------------------------------------
      #### Sum up all rows per id (e.g. sum over xk) ------------------
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
      #### Then sum over the sn splines -------------------------------
      psi_denom <- rowSums(psi_denom)
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      ### by summing across the splines/ columns of psi_t -------------
      w_t <- rowSums(psi_t)
      } else if (errorsY)
      {
      ### P(Y|X,Z)P(Y*|Y,X,Z) -----------------------------------------
      #### Sum up all rows per id (e.g. sum over y) -------------------
      psi_num <- matrix(c(pY_X * pYstar), ncol = 1)
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2))
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      w_t <- psi_t
    }
    # if (max(c(abs(CPPu_t - u_t), abs(CPPw_t - w_t), abs(CPPpsi_t - psi_t))) > 1e-10)
    # {
    #   warning("cpp and R are significantly different!")
    #   browser()
    # }
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    w_t <- .lengthenWT(w_t, n)

    # calculateMu returns exp(-mu) / (1 + exp(-mu))
    muVector <- .calculateMu(theta_design_mat, prev_theta)
    gradient_theta <- .calculateGradient(w_t, n, theta_design_mat, comp_dat_all[, Y], muVector)
    hessian_theta <- .calculateHessian(theta_design_mat, w_t, muVector, n, mus_theta);


    # tic("m-step R")
    # mu <- theta_design_mat %*% prev_theta
    # w_t <- c(rep(1, n), w_t)
    # gradient_theta_R <- matrix(data = c(colSums(w_t * c((comp_dat_all[, Y] - 1 + exp(-mu) / (1 + exp(- mu)))) * theta_design_mat)), ncol = 1)
    # toc()
    #
    #
    # ### ------------------------------------------------------ Gradient
    # ### Hessian -------------------------------------------------------
    # post_multiply <- c((muVector) * (muVector - 1)) * w_t * theta_design_mat
    # hessian_theta_R <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    #
    # if (max(abs(hessian_theta - hessian_theta_R)) > 1e-10)
    # {
    #   warning("cpp and R are significantly different!")
    #   browser()
    # }
    new_theta <- tryCatch(expr = prev_theta - (solve(hessian_theta) %*% gradient_theta),
      error = function(err)
      {
        matrix(NA, nrow = nrow(prev_theta))
        })
    if (any(is.na(new_theta)))
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      # browser()
    }

    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta) < TOL

    ## --------------------------------------------------- Update theta
    ###################################################################
    if (errorsY)
    {
      # w_t is already the proper size

      ## Update gamma using weighted logistic regression ----------------
      muVector <- .calculateMu(gamma_design_mat, prev_gamma)
      gradient_gamma <- .calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[, c(Y_unval)], muVector)
      hessian_gamma <- .calculateHessian(gamma_design_mat, w_t, muVector, n, mus_gamma)

      # mu <- gamma_design_mat %*% prev_gamma
      # gradient_gamma_R <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(- mu) / (1 + exp(- mu)))) * gamma_design_mat)), ncol = 1)

      # if (max(abs(gradient_gamma - gradient_gamma_R)) > 1e-10)
      # {
      #   warning("cpp and R are significantly different! gamma")
      #   browser()
      # }

      # ### ------------------------------------------------------ Gradient
      # ### Hessian -------------------------------------------------------
      # post_multiply <- c(w_t * muVector * (muVector - 1)) * gamma_design_mat
      # hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
      new_gamma <- tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
        error = function(err)
        {
          matrix(NA, nrow = nrow(prev_gamma))
          })
      if (any(is.na(new_gamma)))
      {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      }

      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      ## ---------------- Update gamma using weighted logistic regression
    }
    else
    { gamma_conv <- TRUE }

    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    if (errorsX & errorsY)
    {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num <- p_val_num +
      rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence ---------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
    }
    else if (errorsX)
    {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num <- p_val_num +
      rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence ---------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
    }
    else
    { p_conv <- TRUE }

    ## -------------------------------------------------- Update {p_kj}

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1)
    { CONVERGED <- TRUE }


    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_theta <- new_theta
    if (errorsY)
    { prev_gamma <- new_gamma }

    if (errorsX)
    { prev_p <- new_p }

    #  ------------------------------- Update values for next iteration
  }


  rownames(new_theta) <- c("Intercept", theta_pred)
  if (errorsY)
  {
    rownames(new_gamma) <- c("Intercept", gamma_pred)
  }


  if(!CONVERGED)
  {
    if(it > MAX_ITER)
    {
      CONVERGED_MSG = "MAX_ITER reached"
    }


    return(list(coeff = data.frame(coeff = NA, se = NA),  
      outcome_err_coeff = data.frame(coeff = NA, se = NA),  
      Bspline_coeff = NA, 
      vcov = NA,  
      converged = FALSE,  
      se_converged = NA,  
      converged_msg = CONVERGED_MSG, 
      iterations = it,  
      od_loglik_at_conv = NA))
  }

  if(CONVERGED)
  { CONVERGED_MSG <- "Converged" }

  # ---------------------------------------------- Estimate theta using EM
  if(noSE)
  {
    if (!errorsX)
    {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY)
    {
      new_gamma <- NA
    }

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N,
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

    return(list(coeff = data.frame(coeff = new_theta, se = NA),
      outcome_err_coeff = data.frame(coeff = new_gamma, se = NA),
      Bspline_coeff = cbind(k = comp_dat_val[, "k"], new_p),
      vcov = NA,
      converged = CONVERGED,
      se_converged = NA,
      converged_msg = CONVERGED_MSG,
      iterations = it,
      od_loglik_at_conv = od_loglik_theta))
  }
  else
  {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N <- hn_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    if (!errorsX)
    {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY)
    {
      new_gamma <- NA
    }


    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N,
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

    I_theta <- matrix(od_loglik_theta, nrow = nrow(new_theta), ncol = nrow(new_theta))

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
      I_theta <- matrix(NA, nrow = nrow(new_theta), ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    }
    else
    {
      spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)),
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

      dpt <- matrix(0, nrow = nrow(I_theta), ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta))
      {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
      }

      I_theta <- I_theta + dpt
      # print(I_theta)
    }

    I_theta <- h_N ^ (- 2) * I_theta

    cov_theta <- tryCatch(expr = - solve(I_theta),
      error = function(err)
      {
        matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta))
      }

      )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_theta) < 0)) {
    #   warning("Negative variance estimate. Increase the hn_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED <- FALSE
    # }

    se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
      warning = function(w)
      {
        matrix(NA, nrow = nrow(prev_theta))
        })
    if (any(is.na(se_theta)))
    {
      SE_CONVERGED <- FALSE
    }
    else
    {
      TRUE
    }

    if (verbose)
    {
      message(CONVERGED_MSG)
    }
    
    return(list(coeff = data.frame(coeff = new_theta, se = se_theta),
      outcome_err_coeff = data.frame(coeff = new_gamma, se = NA),
      Bspline_coeff = cbind(k = comp_dat_val[, "k"], new_p),
      vcov = cov_theta,
      converged = CONVERGED,
      se_converged = SE_CONVERGED,
      converged_msg = CONVERGED_MSG,
      iterations = it,
      od_loglik_at_conv = od_loglik_theta))
  }
}

