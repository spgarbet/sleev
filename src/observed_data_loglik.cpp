#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ob_loglik(
    NumericMatrix comp_dat_all,
    int           N,
    int           n,
    IntegerVector theta_pred_cols,
    NumericVector theta,
    int           Y_col,
    IntegerVector gamma_pred_cols,
    NumericVector gamma,
    int           Y_unval_col,
    IntegerVector bspline_cols,
    NumericMatrix p,
    int           k_col
)
{
  int    n_theta    = theta.size();
  int    n_gamma    = gamma.size();
  int    n_bspline  = bspline_cols.size();
  int    m          = p.nrow();
  bool   errorsY    = (n_gamma > 0);
  bool   errorsX    = (n_bspline > 0);
  double loglik     = 0.0;

  // Convert to 0-indexed
  Y_col        -= 1;
  Y_unval_col  -= 1;
  k_col        -= 1;

  // Compute linear predictor and log-likelihood in one pass
  for (int i = 0; i < n; ++i)
  {
    // ===========================================================================
    // Sum over log[P_theta(Yi|Xi)]
    // pY_X <- 1 / (1 + exp(-as.numeric((cbind(int = 1, comp_dat_all[c(1:n), theta_pred]) %*% theta))))
    // pY_X <- ifelse(as.vector(comp_dat_all[c(1:n), c(Y)]) == 0, 1 - pY_X, pY_X)
    // return_loglik <- sum(log(pY_X))

    // Compute linear predictor: intercept + X %*% theta
    double linear_pred = theta[0];  // Intercept
    for (int j = 1; j < n_theta; ++j)
    {
      int pred_col = theta_pred_cols[j - 1] - 1;
      linear_pred += comp_dat_all(i, pred_col) * theta[j];
    }

    // Compute probability: 1 / (1 + exp(-linear_pred))
    double prob = 1.0 / (1.0 + exp(-linear_pred));

    // Adjust probability based on observed Y
    if (comp_dat_all(i, Y_col) == 0.0) prob = 1.0 - prob;

    // Accumulate log-likelihood
    loglik += log(prob);

    // Sum over log[P(Yi*|Xi*,Yi,Xi)]
    // if (errorsY)
    // {
    //   pYstar <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n), gamma_pred]) %*% gamma)))
    //   pYstar <- ifelse(as.vector(comp_dat_all[c(1:n), Y_unval]) == 0, 1 - pYstar, pYstar)
    //   return_loglik <- return_loglik + sum(log(pYstar))
    // }
    if(errorsY)
    {
      double linear_pred_gamma = gamma[0];
      for (int j = 1; j < n_gamma; ++j)
      {
        int pred_col = gamma_pred_cols[j - 1] - 1;
        linear_pred_gamma += comp_dat_all(i, pred_col) * gamma[j];
      }

      double prob_ystar = 1.0 / (1.0 + exp(-linear_pred_gamma));

      if (comp_dat_all(i, Y_unval_col) == 0.0)
        prob_ystar = 1.0 - prob_ystar;

      // Accumulate log-likelihood
      loglik += log(prob_ystar);
    }

    // Sum over I(Xi=xk)Bj(Xi*)log p_kj
    // if (errorsX)
    // {
    //   pX <- p[comp_dat_all[c(1:n), "k"], ]
    //   log_pX <- log(pX)
    //   log_pX[log_pX == -Inf] <- 0
    //   return_loglik <- return_loglik + sum(comp_dat_all[c(1:n), Bspline] * log_pX)
    // }
    if (errorsX)
    {
      // Get k value for this observation (convert from R 1-indexed to C++ 0-indexed)
      int k_val = (int)comp_dat_all(i, k_col) - 1;

      // Sum over Bspline_cols[j] * log(p[k_val, j])
      for (int j = 0; j < n_bspline; ++j)
      {
        double b_val = comp_dat_all(i, bspline_cols[j] - 1);
        double p_val = p(k_val, j);

        // Handle log(0) case: log_pX[log_pX == -Inf] <- 0
        if (p_val > 0.0) loglik += b_val * log(p_val);
      }
    }
  }

  int n_unval = N-n;
  int n_rows_unval = comp_dat_all.nrow() - n;

  // Calculate P_theta(y|x) for all (y,xk)
  // pY_X <- 1 / (1 + exp(-as.numeric((cbind(int = 1, comp_dat_all[-c(1:n), theta_pred]) %*% theta))))
  // pY_X[which(comp_dat_all[-c(1:n), Y] == 0)] <- 1 - pY_X[which(comp_dat_all[-c(1:n), Y] == 0)]
  NumericVector pY_X(n_rows_unval);
  for (int i = 0; i < n_rows_unval; ++i)
  {
    double linear_pred = theta[0];
    for (int j = 1; j < n_theta; ++j)
    {
      int pred_col = theta_pred_cols[j - 1] - 1;
      linear_pred += comp_dat_all(n + i, pred_col) * theta[j];
    }

    double prob = 1.0 / (1.0 + exp(-linear_pred));
    if (comp_dat_all(n + i, Y_col) == 0.0)
      pY_X[i] = 1.0 - prob;
    else
      pY_X[i] = prob;
  }

  // Calculate P(Yi*|Xi*,y,xk) for all (y,xk) if errorsY
  NumericVector pYstar(n_rows_unval);
  if (errorsY)
  {
    for (int i = 0; i < n_rows_unval; ++i)
    {
      double linear_pred_gamma = gamma[0];
      for (int j = 1; j < n_gamma; ++j)
      {
        int pred_col = gamma_pred_cols[j - 1] - 1;
        linear_pred_gamma += comp_dat_all(n + i, pred_col) * gamma[j];
      }

      double prob = 1.0 / (1.0 + exp(-linear_pred_gamma));
      if (comp_dat_all(n + i, Y_unval_col) == 0.0)
        pYstar[i] = 1.0 - prob;
      else
        pYstar[i] = prob;
    }
  }

  // Calculate Bj(Xi*) p_kj for all (k,j) if errorsX
  NumericMatrix pX;
  if (errorsX)
  {
    pX = NumericMatrix(n_rows_unval, n_bspline);
    for (int i = 0; i < n_rows_unval; ++i)
    {
      int k_val = (int)comp_dat_all(n + i, k_col) - 1;
      for (int j = 0; j < n_bspline; ++j)
      {
        pX(i, j) = p(k_val, j);
      }
    }
  }

  // Calculate sum per person and accumulate log-likelihood directly
  if (errorsY && errorsX)
  {
    // Each person has 2*m rows, but they're interleaved
    // Person i has rows: i, i + n_unval, i + 2*n_unval, ...
    int rows_per_person = 2 * m;

    for (int person = 0; person < n_unval; ++person)
    {
      double sum = 0.0;

      for (int j = 0; j < n_bspline; ++j)
      {
        double col_sum = 0.0;

        for (int r = 0; r < rows_per_person; ++r)
        {
          int row_idx = person + r * n_unval;  // Interleaved indexing
          double b_val = comp_dat_all(n + row_idx, bspline_cols[j] - 1);
          double term = pY_X[row_idx] * pYstar[row_idx] * pX(row_idx, j) * b_val;
          col_sum += term;
        }

        sum += col_sum;  // rowSums: sum across columns
      }

      // Accumulate log directly (handle -Inf case)
      if (sum > 0.0) loglik += log(sum);
    }
  }
  else if (errorsY)
  {
    // Each person has 2 rows, interleaved
    int rows_per_person = 2;

    for (int person = 0; person < n_unval; ++person)
    {
      double sum = 0.0;

      for (int r = 0; r < rows_per_person; ++r)
      {
        int row_idx = person + r * n_unval;  // Interleaved indexing
        sum += pY_X[row_idx] * pYstar[row_idx];
      }

      // Accumulate log directly (handle -Inf case)
      if (sum > 0.0) loglik += log(sum);
    }
  }
  else if (errorsX)
  {
    // Each person has m rows, interleaved
    int rows_per_person = m;

    for (int person = 0; person < n_unval; ++person)
    {
      double sum = 0.0;

      for (int j = 0; j < n_bspline; ++j)
      {
        double col_sum = 0.0;

        for (int r = 0; r < rows_per_person; ++r)
        {
          int row_idx = person + r * n_unval;  // Interleaved indexing
          double b_val = comp_dat_all(n + row_idx, bspline_cols[j] - 1);
          double term = pY_X[row_idx] * pX(row_idx, j) * b_val;
          col_sum += term;
        }

        sum += col_sum;  // rowSums: sum across columns
      }

      // Accumulate log directly (handle -Inf case)
      if (sum > 0.0) loglik += log(sum);
    }
  }

  return loglik;
}
