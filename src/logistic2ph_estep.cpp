#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List logistic2ph_estep_cpp(
    const arma::mat& theta_design,
    const arma::mat& gamma_design,
    const arma::mat& prev_theta,
    const arma::mat& prev_gamma,
    const arma::mat& prev_p,
    const arma::mat& Bspline,
    const arma::colvec& Y,
    const arma::colvec& Y_unval,
    int n,
    int N,
    int m,
    bool errorsY
)
{
  const int Nu = N - n;
  const int n_spline = Bspline.n_cols;
  const int mult = errorsY ? 2 : 1;
  const int total_rows = mult * m * Nu;

  // Validate dimensions
  if (Bspline.n_rows != total_rows) {
    Rcpp::stop("Dimension mismatch: Bspline has %d rows, expected %d",
               Bspline.n_rows, total_rows);
  }

  // p(Y | X) - compute for ALL rows in comp_dat_unval
  arma::colvec mu_theta = theta_design.rows(n, theta_design.n_rows - 1) * prev_theta;
  arma::colvec pY_X = 1.0 / (1.0 + arma::exp(-mu_theta));
  for (int i = 0; i < total_rows; ++i) {
    if (Y[i] == 0) pY_X[i] = 1.0 - pY_X[i];
  }

  // p(Y* | X*, Y, X)
  arma::colvec pYstar(total_rows, arma::fill::ones);
  if (errorsY)
  {
    arma::colvec mu_gamma = gamma_design.rows(n, gamma_design.n_rows - 1) * prev_gamma;
    for (int i = 0; i < total_rows; ++i) {
      double p = 1.0 / (1.0 + std::exp(-mu_gamma[i]));
      pYstar[i] = (Y_unval[i] == 0) ? (1.0 - p) : p;
    }
  }

  // Build p(X|X*) and psi_num
  arma::mat psi_num(total_rows, n_spline);
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    for (int k = 0; k < m; ++k) {
      for (int i = 0; i < Nu; ++i) {
        int row_idx = y_idx * (m * Nu) + k * Nu + i;
        double py = pY_X[row_idx] * pYstar[row_idx];

        for (int j = 0; j < n_spline; ++j) {
          psi_num(row_idx, j) = py * prev_p(k, j) * Bspline(row_idx, j);
        }
      }
    }
  }

  // Compute denominator by subject
  arma::colvec psi_denom(Nu, arma::fill::zeros);
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    for (int k = 0; k < m; ++k) {
      for (int i = 0; i < Nu; ++i) {
        int row_idx = y_idx * (m * Nu) + k * Nu + i;
        psi_denom[i] += arma::accu(psi_num.row(row_idx));
      }
    }
  }

  for (int i = 0; i < Nu; ++i) {
    if (psi_denom[i] == 0.0) psi_denom[i] = 1.0;
  }

  // Compute psi_t
  arma::mat psi_t(total_rows, n_spline);
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    for (int k = 0; k < m; ++k) {
      for (int i = 0; i < Nu; ++i) {
        int row_idx = y_idx * (m * Nu) + k * Nu + i;
        psi_t.row(row_idx) = psi_num.row(row_idx) / psi_denom[i];
      }
    }
  }

  // Compute w_t
  arma::colvec w_t = arma::sum(psi_t, 1);

  // Compute u_t
  arma::mat u_t = psi_t.rows(0, m * Nu - 1);
  if (errorsY) {
    u_t += psi_t.rows(m * Nu, total_rows - 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("w_t") = w_t,
    Rcpp::Named("u_t") = u_t
  );
}
