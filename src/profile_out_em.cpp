// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "sleev_api.h"

// Complete EM algorithm for profiling out nuisance parameters
// [[Rcpp::export]]
Rcpp::List profile_out_em(
    const  arma::mat&  theta_design_mat,
    const  arma::vec&  theta,
    const  arma::mat&  comp_dat_all,
    int                Y_col,
    const  arma::mat&  gamma_design_mat,
    const  arma::vec&  gamma0,
    int                Y_unval_col,
    const  arma::mat&  p0,
    const  arma::uvec& Bspline_cols,
    const  arma::mat&  p_val_num,
    int                n,
    int                N,
    int                m,
    bool               errorsY,
    double             TOL,
    int                MAX_ITER,
    bool               verbose)
{
  // // Access the exported functions from sleev package
  // Environment sleev("package:sleev");
  // Function    pYstarCalc        = sleev[".pYstarCalc"];
  // Function    lengthenWT        = sleev[".lengthenWT"];
  // Function    calculateMu       = sleev[".calculateMu"];
  // Function    calculateGradient = sleev[".calculateGradient"];
  // Function    calculateHessian  = sleev[".calculateHessian"];

  int         N_minus_n         = N - n;
  arma::vec   prev_gamma        = gamma0;
  arma::mat   prev_p            = p0;
  bool        CONVERGED         = false;
  std::string CONVERGED_MSG     = "Unknown";
  int         it                = 1;

  // Pre-allocate for M-step
  arma::vec mus_gamma;
  if (errorsY)
  {
    mus_gamma = arma::vec(gamma_design_mat.n_rows * prev_gamma.n_elem);
  }

  // Variables to hold results
  arma::mat psi_t;
  arma::vec w_t;
  arma::mat u_t;
  arma::vec new_gamma;
  arma::mat new_p;

  Rcpp::Rcout << "Beginning Profile EM Loop" << std::endl;

  while (it <= MAX_ITER && !CONVERGED)
  {
    if (verbose && (it % 10 == 0 || it < 5))
    {
      Rcpp::Rcout << "Iteration: " << (it + 1) << std::endl;
    }

     ///////////////////////////////////////////////////////////////////
    // E Step

    arma::vec pY_X = pYstarCalc(theta_design_mat, 1, theta,
                                comp_dat_all, Y_col);

    // Calculate pYstar: P(Y*|X*,Y,X)
    arma::vec pYstar;
    if (errorsY)
    {
      pYstar = pYstarCalc(gamma_design_mat, n + 1, prev_gamma,
                          comp_dat_all, Y_unval_col);
    }
    else
    {
      pYstar = arma::vec(pY_X.n_elem, arma::fill::ones);
    }

    // Calculate pX: P(X|X*)
    // Extract B-spline columns
    int n_bspline_cols = Bspline_cols.n_elem;
    arma::mat bspline_data(comp_dat_all.n_rows - n, n_bspline_cols);
    for (int i = 0; i < bspline_data.n_rows; ++i)
    {
      for (int j = 0; j < n_bspline_cols; ++j)
      {
        bspline_data(i, j) = comp_dat_all(n + i, Bspline_cols[j]);
      }
    }

    // Create pX matrix
    int pX_rows;
    if (errorsY)
    {
      pX_rows = m * N_minus_n * 2;
    }
    else
    {
      pX_rows = m * N_minus_n;
    }

    arma::mat pX(pX_rows, n_bspline_cols);

    // Fill pX with reordered prev_p multiplied by B-spline data
    if (errorsY)
    {
      // rep(rep(seq(1, m), each = (N - n)), times = 2)
      for (int rep_idx = 0; rep_idx < 2; ++rep_idx)
      {
        for (int k = 0; k < m; ++k)
        {
          for (int i = 0; i < N_minus_n; ++i)
          {
            int pX_row = rep_idx * m * N_minus_n + k * N_minus_n + i;
            for (int j = 0; j < n_bspline_cols; ++j)
            {
              pX(pX_row, j) = prev_p(k, j) * bspline_data(i, j);
            }
          }
        }
      }
    }
    else
    {
      // rep(seq(1, m), each = (N - n))
      for (int k = 0; k < m; ++k)
      {
        for (int i = 0; i < N_minus_n; ++i)
        {
          int pX_row = k * N_minus_n + i;
          for (int j = 0; j < n_bspline_cols; ++j)
          {
            pX(pX_row, j) = prev_p(k, j) * bspline_data(i, j);
          }
        }
      }
    }

    // Calculate conditional expectations
    int n_rows = pX.n_rows;
    int n_cols = pX.n_cols;

    // Allocate matrices
    arma::mat psi_num(n_rows, n_cols);
    psi_t = arma::mat(n_rows, n_cols);
    w_t   = arma::vec(n_rows);
    u_t   = arma::mat(m * N_minus_n, n_cols);

    // Calculate psi_num: c(pY_X * pYstar) * pX
    int pYstar_len = pYstar.n_elem;

    for (int i = 0; i < n_rows; ++i)
    {
      // Handle R-style recycling for pYstar
      int pYstar_idx = (pYstar_len == 1) ? 0 : i;
      double weight = pY_X[i] * pYstar[pYstar_idx];

      for (int j = 0; j < n_cols; ++j)
      {
        psi_num(i, j) = weight * pX(i, j);
      }
    }

    // Calculate psi_denom by summing over groups using rowsum logic
    arma::vec psi_denom(N_minus_n, arma::fill::zeros);

    for (int row = 0; row < n_rows; ++row)
    {
      int group_idx = row % N_minus_n;

      for (int j = 0; j < n_cols; ++j)
      {
        psi_denom[group_idx] += psi_num(row, j);
      }
    }

    // Avoid division by zero
    for (int i = 0; i < N_minus_n; ++i)
    {
      if (psi_denom[i] == 0.0)
      {
        psi_denom[i] = 1.0;
      }
    }

    // Calculate psi_t: psi_num / psi_denom
    for (int row = 0; row < n_rows; ++row)
    {
      int group_idx = row % N_minus_n;
      for (int j = 0; j < n_cols; ++j)
      {
        psi_t(row, j) = psi_num(row, j) / psi_denom[group_idx];
      }
    }

    // Calculate w_t by summing across columns (rowSums)
    for (int i = 0; i < n_rows; ++i)
    {
      w_t[i] = 0.0;
      for (int j = 0; j < n_cols; ++j)
      {
        w_t[i] += psi_t(i, j);
      }
    }

    // Calculate u_t
    int u_rows = m * N_minus_n;

    for (int i = 0; i < u_rows; ++i)
    {
      for (int j = 0; j < n_cols; ++j)
      {
        u_t(i, j) = psi_t(i, j);
        if (errorsY)
        {
          u_t(i, j) += psi_t(i + u_rows, j);
        }
      }
    }

     ///////////////////////////////////////////////////////////////////
    // M Step

    // Update gamma using weighted logistic regression
    arma::vec gamma_conv;
    if (errorsY)
    {
      w_t = lengthenWT(w_t, n, true);
      arma::vec muVector = calculateMu(gamma_design_mat, prev_gamma);

      // Extract Y_unval column
      arma::vec Y_unval_vec(comp_dat_all.n_rows);
      for (int i = 0; i < comp_dat_all.n_rows; ++i)
      {
        Y_unval_vec[i] = comp_dat_all(i, Y_unval_col);
      }

      arma::vec gradient_gamma = calculateGradient(w_t, n, gamma_design_mat,
                                                   Y_unval_vec, muVector, false);
      arma::mat hessian_gamma  = calculateHessian(gamma_design_mat, w_t,
                                                  muVector, n, mus_gamma, false);

      // Try to solve: new_gamma = prev_gamma - solve(hessian_gamma) %*% gradient_gamma
      try
      {
        arma::vec update = arma::solve(hessian_gamma, gradient_gamma);
        new_gamma        = prev_gamma - update;

        // Check if any values are NA/inf
        bool has_na = false;
        for (int i = 0; i < new_gamma.n_elem; ++i)
        {
          if (!arma::is_finite(new_gamma[i]))
          {
            has_na = true;
            break;
          }
        }

        if (has_na)
        {
          // Create formula dynamically would be complex, so we throw an error
          // In practice, this fallback path is rarely hit
          stop("Matrix inversion failed and GLM fallback not implemented in C++");
        }
      }
      catch (...)
      {
        stop("Matrix inversion failed in M-step");
      }

      // Check for convergence
      gamma_conv = arma::abs(new_gamma - prev_gamma);
    }
    else
    {
      new_gamma  = arma::vec();  // Empty vector
      gamma_conv = arma::vec(1, arma::fill::zeros);  // Will pass convergence test
    }

    // Update {p_kj}
    // new_p_num <- p_val_num + rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    arma::mat new_p_num = p_val_num;

    // Add u_t grouped by rep(seq(1, m), each = (N - n))
    for (int i = 0; i < u_t.n_rows; ++i)
    {
      int group_idx = i / N_minus_n;  // Integer division gives the group
      for (int j = 0; j < u_t.n_cols; ++j)
      {
        new_p_num(group_idx, j) += u_t(i, j);
      }
    }

    // new_p <- t(t(new_p_num) / colSums(new_p_num))
    // First calculate colSums
    arma::rowvec col_sums = arma::sum(new_p_num, 0);

    // Divide each column by its sum
    new_p = arma::mat(new_p_num.n_rows, new_p_num.n_cols);
    for (int j = 0; j < new_p.n_cols; ++j)
    {
      new_p.col(j) = new_p_num.col(j) / col_sums[j];
    }

    // Check for convergence
    arma::vec p_conv   = arma::abs(arma::vectorise(new_p - prev_p));

    // Check overall convergence
    arma::vec all_conv = arma::join_cols(gamma_conv, p_conv);

    bool all_true = arma::all(all_conv < TOL);

    if (all_true)
    {
      CONVERGED = true;
    }

    // Update for next iteration
    ++it;
    prev_gamma = new_gamma;
    prev_p = new_p;
  }

  if (it == MAX_ITER + 1 && !CONVERGED)
  {
    CONVERGED_MSG = "MAX_ITER reached";
    if (errorsY)
    {
      new_gamma = arma::vec(gamma0.n_elem);
      new_gamma.fill(arma::datum::nan);
    }
    else
    {
      new_gamma = arma::vec();  // Empty vector
    }
    new_p = arma::mat(p0.n_rows, p0.n_cols);
    new_p.fill(arma::datum::nan);
  }

  if (CONVERGED)
  {
    CONVERGED_MSG = "converged";
  }

  Rcpp::Rcout << "Ending Profile EM Loop" << std::endl;

  return List::create(
    Named("psi_at_conv")   = psi_t,
    Named("gamma_at_conv") = new_gamma,
    Named("p_at_conv")     = new_p,
    Named("converged")     = CONVERGED,
    Named("converged_msg") = CONVERGED_MSG
  );
}
