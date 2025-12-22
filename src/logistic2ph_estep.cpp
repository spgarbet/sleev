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

  // ===================================================================
  // Allocate final output matrices ONCE - reuse throughout
  // ===================================================================
  arma::mat psi_t(total_rows, n_spline);  // Will be our main workspace
  arma::colvec psi_denom(Nu, arma::fill::zeros);

  // ===================================================================
  // STEP 1: Compute probabilities using raw pointer to first column
  // ===================================================================
  // Use direct pointer access to first column as workspace
  double* py_ptr = psi_t.colptr(0);

  // Compute p(Y | X) directly in workspace
  arma::colvec mu_theta = theta_design.rows(n, theta_design.n_rows - 1) * prev_theta;
  for (int i = 0; i < total_rows; ++i) {
    py_ptr[i] = 1.0 / (1.0 + std::exp(-mu_theta[i]));
    if (Y[i] == 0) py_ptr[i] = 1.0 - py_ptr[i];
  }

  // Multiply by p(Y*|X*,Y,X) if needed (in-place)
  if (errorsY) {
    arma::colvec mu_gamma = gamma_design.rows(n, gamma_design.n_rows - 1) * prev_gamma;
    for (int i = 0; i < total_rows; ++i) {
      double p = 1.0 / (1.0 + std::exp(-mu_gamma[i]));
      py_ptr[i] *= (Y_unval[i] == 0) ? (1.0 - p) : p;
    }
  }

  // ===================================================================
  // STEP 2: Compute psi_num directly in psi_t (no intermediate matrix)
  // ===================================================================
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    int y_offset = y_idx * m * Nu;

    for (int k = 0; k < m; ++k) {
      int k_offset = y_offset + k * Nu;

      for (int i = 0; i < Nu; ++i) {
        int row_idx = k_offset + i;
        double py_val = py_ptr[row_idx];

        // Compute psi_num for all splines j in this row
        // Direct pointer access for better cache performance
        for (int j = 0; j < n_spline; ++j) {
          psi_t(row_idx, j) = py_val * prev_p(k, j) * Bspline(row_idx, j);
        }
      }
    }
  }

  // ===================================================================
  // STEP 3: Compute denominator (single pass, no extra allocations)
  // ===================================================================
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    for (int k = 0; k < m; ++k) {
      int k_offset = y_idx * (m * Nu) + k * Nu;

      for (int i = 0; i < Nu; ++i) {
        int row_idx = k_offset + i;

        // Sum across splines for this row
        double row_sum = 0.0;
        for (int j = 0; j < n_spline; ++j) {
          row_sum += psi_t(row_idx, j);
        }

        psi_denom[i] += row_sum;
      }
    }
  }

  // Avoid division by zero
  for (int i = 0; i < Nu; ++i) {
    if (psi_denom[i] == 0.0) psi_denom[i] = 1.0;
  }

  // ===================================================================
  // STEP 4: Normalize psi_t IN-PLACE (no new allocation)
  // ===================================================================
  for (int y_idx = 0; y_idx < mult; ++y_idx) {
    for (int k = 0; k < m; ++k) {
      int k_offset = y_idx * (m * Nu) + k * Nu;

      for (int i = 0; i < Nu; ++i) {
        int row_idx = k_offset + i;
        double denom = psi_denom[i];

        // Divide entire row by denominator
        for (int j = 0; j < n_spline; ++j) {
          psi_t(row_idx, j) /= denom;
        }
      }
    }
  }

  // ===================================================================
  // STEP 5: Compute w_t (sum across columns)
  // ===================================================================
  arma::colvec w_t(total_rows, arma::fill::zeros);

  for (int i = 0; i < total_rows; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n_spline; ++j) {
      sum += psi_t(i, j);
    }
    w_t[i] = sum;
  }

  // ===================================================================
  // STEP 6: Compute u_t (collapse Y dimension)
  // ===================================================================
  arma::mat u_t(m * Nu, n_spline);

  if (errorsY) {
    // Add Y=0 and Y=1 contributions
    for (int k = 0; k < m; ++k) {
      for (int i = 0; i < Nu; ++i) {
        int out_row = k * Nu + i;
        int row_y0 = k * Nu + i;
        int row_y1 = (m * Nu) + k * Nu + i;

        for (int j = 0; j < n_spline; ++j) {
          u_t(out_row, j) = psi_t(row_y0, j) + psi_t(row_y1, j);
        }
      }
    }
  } else {
    // Just copy first m*Nu rows
    u_t = psi_t.rows(0, m * Nu - 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("w_t") = w_t,
    Rcpp::Named("u_t") = u_t
  );
}
