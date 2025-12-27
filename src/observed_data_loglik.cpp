#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double compute_validated_y_loglik(
    NumericMatrix comp_dat_all,
    int           n,
    IntegerVector theta_pred_cols,
    NumericVector theta,
    int           Y_col,
    IntegerVector gamma_pred_cols,
    NumericVector gamma,
    int           Y_unval_col
)
{
  int    n_theta    = theta.size();
  int    n_gamma    = gamma.size();
  bool   errorsY    = (n_gamma > 0);
  double loglik     = 0.0;

  // Convert to 0-indexed
  Y_col       -= 1;
  Y_unval_col -= 1;

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
  }

  return loglik;
}

// [[Rcpp::export]]
double observed_data_loglik_cpp(
    int N,
    int n,
    IntegerVector Y_unval,
    CharacterVector Y,
    IntegerVector X_unval,
    CharacterVector X,
    CharacterVector Z,
    CharacterVector Bspline,
    NumericMatrix comp_dat_all,
    CharacterVector theta_pred,
    CharacterVector gamma_pred,
    NumericVector theta,
    NumericVector gamma,
    NumericMatrix p
)
{
  // Error settings based on arguments provided
  bool errorsY = (Y_unval.size() > 0);
  bool errorsX = (X_unval.size() > 0);

  int m = p.nrow();
  double return_loglik = 0.0;

  // Get column names once
  CharacterVector col_names = Rcpp::colnames(comp_dat_all);

  // Helper function to find column index by name
  auto find_col = [&col_names](const String& name) -> int
  {
    for (int i = 0; i < col_names.size(); i++)
    {
      if (col_names[i] == name)
      {
        return i;
      }
    }
    return -1;
  };

  // Find Y column index
  int Y_col = -1;
  if (Y.size() > 0)
  {
    Y_col = find_col(Y[0]);
  }

  // ============================================================================
  // FOR VALIDATED SUBJECTS (rows 0 to n-1 in C++ indexing)
  // ============================================================================

  // Sum over log[P_theta(Yi|Xi)]
  // Build design matrix: cbind(int=1, comp_dat_all[1:n, theta_pred])
  int n_theta = theta.size();
  NumericMatrix X_theta_val(n, n_theta);

  // Intercept
  for (int i = 0; i < n; i++)
  {
    X_theta_val(i, 0) = 1.0;
  }

  // Predictors
  for (int j = 1; j < n_theta; j++)
  {
    int col_idx = find_col(theta_pred[j - 1]);
    for (int i = 0; i < n; i++)
    {
      X_theta_val(i, j) = comp_dat_all(i, col_idx);
    }
  }

  // Compute P(Y|X) for validated
  for (int i = 0; i < n; i++)
  {
    double linear_pred = 0.0;
    for (int j = 0; j < n_theta; j++)
    {
      linear_pred += X_theta_val(i, j) * theta[j];
    }

    double pY_X = 1.0 / (1.0 + exp(-linear_pred));
    if (comp_dat_all(i, Y_col) == 0.0)
    {
      pY_X = 1.0 - pY_X;
    }

    return_loglik += log(pY_X);
  }

  // Sum over log[P(Yi*|Xi*,Yi,Xi)]
  if (errorsY)
  {
    int Y_unval_col = Y_unval[0] - 1; // Convert from R 1-indexed to C++ 0-indexed
    int n_gamma = gamma.size();
    NumericMatrix X_gamma_val(n, n_gamma);

    // Intercept
    for (int i = 0; i < n; i++)
    {
      X_gamma_val(i, 0) = 1.0;
    }

    // Predictors
    for (int j = 1; j < n_gamma; j++)
    {
      int col_idx = find_col(gamma_pred[j - 1]);
      for (int i = 0; i < n; i++)
      {
        X_gamma_val(i, j) = comp_dat_all(i, col_idx);
      }
    }

    // Compute P(Y*|X*,Y,X) for validated
    for (int i = 0; i < n; i++)
    {
      double linear_pred = 0.0;
      for (int j = 0; j < n_gamma; j++)
      {
        linear_pred += X_gamma_val(i, j) * gamma[j];
      }

      double pYstar = 1.0 / (1.0 + exp(-linear_pred));
      if (comp_dat_all(i, Y_unval_col) == 0.0)
      {
        pYstar = 1.0 - pYstar;
      }

      return_loglik += log(pYstar);
    }
  }

  // Sum over I(Xi=xk)Bj(Xi*)log p_kj
  if (errorsX)
  {
    int k_col = find_col("k");
    int n_bspline = Bspline.size();

    // Find Bspline column indices
    IntegerVector bspline_cols(n_bspline);
    for (int j = 0; j < n_bspline; j++)
    {
      bspline_cols[j] = find_col(Bspline[j]);
    }

    for (int i = 0; i < n; i++)
    {
      int k_val = (int)comp_dat_all(i, k_col) - 1; // Convert to 0-indexed

      for (int j = 0; j < n_bspline; j++)
      {
        double b_val = comp_dat_all(i, bspline_cols[j]);
        double p_val = p(k_val, j);

        if (p_val > 0.0)
        {
          return_loglik += b_val * log(p_val);
        }
      }
    }
  }

  // ============================================================================
  // FOR UNVALIDATED SUBJECTS (rows n to N-1 in C++ indexing)
  // ============================================================================

  int n_unval = N - n;

  // Calculate P_theta(y|x) for all (y,xk)
  int n_rows_unval = comp_dat_all.nrow() - n;
  NumericVector pY_X_unval(n_rows_unval);

  // Build design matrix for unvalidated
  NumericMatrix X_theta_unval(n_rows_unval, n_theta);
  for (int i = 0; i < n_rows_unval; i++)
  {
    X_theta_unval(i, 0) = 1.0;
    for (int j = 1; j < n_theta; j++)
    {
      int col_idx = find_col(theta_pred[j - 1]);
      X_theta_unval(i, j) = comp_dat_all(n + i, col_idx);
    }
  }

  // Compute P(Y|X)
  for (int i = 0; i < n_rows_unval; i++)
  {
    double linear_pred = 0.0;
    for (int j = 0; j < n_theta; j++)
    {
      linear_pred += X_theta_unval(i, j) * theta[j];
    }

    double prob = 1.0 / (1.0 + exp(-linear_pred));
    if (comp_dat_all(n + i, Y_col) == 0.0)
    {
      pY_X_unval[i] = 1.0 - prob;
    }
    else
    {
      pY_X_unval[i] = prob;
    }
  }

  // Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  NumericVector pYstar_unval(n_rows_unval);
  if (errorsY)
  {
    int Y_unval_col = Y_unval[0] - 1;
    int n_gamma = gamma.size();
    NumericMatrix X_gamma_unval(n_rows_unval, n_gamma);

    for (int i = 0; i < n_rows_unval; i++)
    {
      X_gamma_unval(i, 0) = 1.0;
      for (int j = 1; j < n_gamma; j++)
      {
        int col_idx = find_col(gamma_pred[j - 1]);
        X_gamma_unval(i, j) = comp_dat_all(n + i, col_idx);
      }
    }

    for (int i = 0; i < n_rows_unval; i++)
    {
      double linear_pred = 0.0;
      for (int j = 0; j < n_gamma; j++)
      {
        linear_pred += X_gamma_unval(i, j) * gamma[j];
      }

      double prob = 1.0 / (1.0 + exp(-linear_pred));
      if (comp_dat_all(n + i, Y_unval_col) == 0.0)
      {
        pYstar_unval[i] = 1.0 - prob;
      }
      else
      {
        pYstar_unval[i] = prob;
      }
    }
  }

  // Calculate Bj(Xi*) p_kj for all (k,j)
  NumericMatrix pX_unval;
  if (errorsX)
  {
    int k_col = find_col("k");
    pX_unval = NumericMatrix(n_rows_unval, p.ncol());

    for (int i = 0; i < n_rows_unval; i++)
    {
      int k_val = (int)comp_dat_all(n + i, k_col) - 1; // Convert to 0-indexed
      for (int j = 0; j < p.ncol(); j++)
      {
        pX_unval(i, j) = p(k_val, j);
      }
    }
  }

  // Calculate sum of P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj
  // Need to sum over rows that belong to same person
  NumericVector person_sum(n_unval);

  if (errorsY && errorsX)
  {
    int n_bspline = Bspline.size();
    IntegerVector bspline_cols(n_bspline);
    for (int j = 0; j < n_bspline; j++)
    {
      bspline_cols[j] = find_col(Bspline[j]);
    }

    // Each person has 2*m rows in the expanded data
    int rows_per_person = 2 * m;

    for (int person = 0; person < n_unval; person++)
    {
      double sum = 0.0;
      int start_row = person * rows_per_person;

      for (int r = 0; r < rows_per_person; r++)
      {
        int row_idx = start_row + r;

        // For each Bspline column j:
        // sum over: pY_X * pYstar * pX[,j] * Bspline[,j]
        for (int j = 0; j < n_bspline; j++)
        {
          double b_val = comp_dat_all(n + row_idx, bspline_cols[j]);
          double p_val = pX_unval(row_idx, j);
          double term = pY_X_unval[row_idx] * pYstar_unval[row_idx] * p_val * b_val;
          sum += term;
        }
      }

      person_sum[person] = sum;
    }
  }
  else if (errorsY)
  {
    // Each person has 2 rows in the expanded data
    int rows_per_person = 2;

    for (int person = 0; person < n_unval; person++)
    {
      double sum = 0.0;
      int start_row = person * rows_per_person;

      for (int r = 0; r < rows_per_person; r++)
      {
        int row_idx = start_row + r;
        sum += pY_X_unval[row_idx] * pYstar_unval[row_idx];
      }

      person_sum[person] = sum;
    }
  }
  else if (errorsX)
  {
    int n_bspline = Bspline.size();
    IntegerVector bspline_cols(n_bspline);
    for (int j = 0; j < n_bspline; j++)
    {
      bspline_cols[j] = find_col(Bspline[j]);
    }

    // Each person has m rows in the expanded data
    int rows_per_person = m;

    for (int person = 0; person < n_unval; person++)
    {
      double sum = 0.0;
      int start_row = person * rows_per_person;

      for (int r = 0; r < rows_per_person; r++)
      {
        int row_idx = start_row + r;

        // For each Bspline column j:
        // sum over: pY_X * pX[,j] * Bspline[,j]
        for (int j = 0; j < n_bspline; j++)
        {
          double b_val = comp_dat_all(n + row_idx, bspline_cols[j]);
          double p_val = pX_unval(row_idx, j);
          double term = pY_X_unval[row_idx] * p_val * b_val;
          sum += term;
        }
      }

      person_sum[person] = sum;
    }
  }

  // Sum log of person sums
  for (int person = 0; person < n_unval; person++)
  {
    if (person_sum[person] > 0.0)
    {
      return_loglik += log(person_sum[person]);
    }
  }

  return return_loglik;
}
