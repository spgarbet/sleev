#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function: compute mu = exp(-X*beta) / (1 + exp(-X*beta))
// This is used in the gradient and Hessian calculations
inline void compute_mu(
    const arma::mat& design_mat,
    const arma::colvec& beta,
    arma::colvec& mu_out)
{
  // mu = exp(-X*beta) / (1 + exp(-X*beta))
  arma::colvec xb = design_mat * beta;
  for (size_t i = 0; i < xb.n_elem; ++i) {
    double exp_neg_xb = std::exp(-xb[i]);
    mu_out[i] = exp_neg_xb / (1.0 + exp_neg_xb);
  }
}

// Compute gradient: sum_i w_i * (y_i - (1 - mu_i)) * x_i
// = sum_i w_i * (y_i - 1 + mu_i) * x_i
inline void compute_gradient(
    const arma::mat& design_mat,
    const arma::colvec& w,
    const arma::colvec& y,
    const arma::colvec& mu,
    int n_validated,
    arma::colvec& gradient_out)
{
  int n_total = design_mat.n_rows;
  int n_params = design_mat.n_cols;

  gradient_out.zeros(n_params);

  // For validated subjects (first n rows), weight is 1
  for (int i = 0; i < n_validated; ++i) {
    double residual = y[i] - 1.0 + mu[i];
    for (int j = 0; j < n_params; ++j) {
      gradient_out[j] += residual * design_mat(i, j);
    }
  }

  // For unvalidated subjects, use weight w[i]
  for (int i = n_validated; i < n_total; ++i) {
    double residual = y[i] - 1.0 + mu[i];
    double weight = w[i];
    for (int j = 0; j < n_params; ++j) {
      gradient_out[j] += weight * residual * design_mat(i, j);
    }
  }
}

// Compute Hessian: sum_i w_i * mu_i * (1 - mu_i) * x_i * x_i^T
inline void compute_hessian(
    const arma::mat& design_mat,
    const arma::colvec& w,
    const arma::colvec& mu,
    int n_validated,
    arma::mat& hessian_out)
{
  int n_total = design_mat.n_rows;
  int n_params = design_mat.n_cols;

  hessian_out.zeros(n_params, n_params);

  // For validated subjects (first n rows), weight is 1
  for (int i = 0; i < n_validated; ++i) {
    double mu_i = mu[i];
    double scale = mu_i * (1.0 - mu_i);

    for (int j = 0; j < n_params; ++j) {
      double x_ij = design_mat(i, j);
      for (int k = j; k < n_params; ++k) {
        hessian_out(j, k) += scale * x_ij * design_mat(i, k);
      }
    }
  }

  // For unvalidated subjects, use weight w[i]
  for (int i = n_validated; i < n_total; ++i) {
    double mu_i = mu[i];
    double weight = w[i];
    double scale = weight * mu_i * (1.0 - mu_i);

    for (int j = 0; j < n_params; ++j) {
      double x_ij = design_mat(i, j);
      for (int k = j; k < n_params; ++k) {
        hessian_out(j, k) += scale * x_ij * design_mat(i, k);
      }
    }
  }

  // Fill lower triangle (Hessian is symmetric)
  for (int j = 0; j < n_params; ++j) {
    for (int k = j + 1; k < n_params; ++k) {
      hessian_out(k, j) = hessian_out(j, k);
    }
  }
}

// Newton-Raphson update: beta_new = beta_old - H^{-1} * g
// Returns true if successful, false if Hessian is singular
inline bool newton_raphson_step(
    const arma::colvec& beta_old,
    const arma::colvec& gradient,
    const arma::mat& hessian,
    arma::colvec& beta_new)
{
  // Try to solve H * delta = -gradient
  arma::colvec delta;
  bool success = arma::solve(delta, hessian, -gradient, arma::solve_opts::likely_sympd);

  if (success) {
    beta_new = beta_old + delta;
    return true;
  } else {
    beta_new = beta_old;
    return false;
  }
}

// Update p (B-spline coefficients)
// new_p_num = p_val_num + rowsum(u_t grouped by k)
// new_p = new_p_num / colSums(new_p_num)
inline void update_p(
    const arma::mat& p_val_num,
    const arma::mat& u_t,
    int m,
    int Nu,
    arma::mat& p_new)
{
  int n_spline = p_val_num.n_cols;

  // Compute new_p_num = p_val_num + rowsum(u_t, group by k)
  arma::mat new_p_num = p_val_num;

  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < Nu; ++i) {
      int row_idx = k * Nu + i;
      for (int j = 0; j < n_spline; ++j) {
        new_p_num(k, j) += u_t(row_idx, j);
      }
    }
  }

  // Normalize: p_new = new_p_num / colSums(new_p_num)
  arma::rowvec col_sums = arma::sum(new_p_num, 0);

  for (int k = 0; k < m; ++k) {
    for (int j = 0; j < n_spline; ++j) {
      p_new(k, j) = new_p_num(k, j) / col_sums[j];
    }
  }
}

// Complete M-step update
// [[Rcpp::export]]
Rcpp::List logistic2ph_mstep(
    const arma::mat& theta_design_mat,
    const arma::mat& gamma_design_mat,
    const arma::colvec& Y_all,
    const arma::colvec& Y_unval_all,
    const arma::colvec& w_t,
    const arma::mat& u_t,
    const arma::mat& p_val_num,
    const arma::colvec& prev_theta,
    const arma::colvec& prev_gamma,
    const arma::mat& prev_p,
    int n,
    int N,
    int m,
    bool errorsY,
    double TOL
)
{
  int Nu = N - n;
  int n_params_theta = theta_design_mat.n_cols;
  int n_params_gamma = gamma_design_mat.n_cols;

  // ===================================================================
  // Update theta
  // ===================================================================
  arma::colvec mu_theta(theta_design_mat.n_rows);
  arma::colvec gradient_theta(n_params_theta);
  arma::mat hessian_theta(n_params_theta, n_params_theta);
  arma::colvec new_theta(n_params_theta);

  compute_mu(theta_design_mat, prev_theta, mu_theta);
  compute_gradient(theta_design_mat, w_t, Y_all, mu_theta, n, gradient_theta);
  compute_hessian(theta_design_mat, w_t, mu_theta, n, hessian_theta);

  bool theta_success = newton_raphson_step(prev_theta, gradient_theta, hessian_theta, new_theta);

  // Check convergence
  arma::colvec theta_diff = arma::abs(new_theta - prev_theta);
  bool theta_conv = arma::all(theta_diff < TOL);

  // ===================================================================
  // Update gamma (if errorsY)
  // ===================================================================
  arma::colvec new_gamma;
  bool gamma_success = true;
  bool gamma_conv = true;

  if (errorsY)
  {
    arma::colvec mu_gamma(gamma_design_mat.n_rows);
    arma::colvec gradient_gamma(n_params_gamma);
    arma::mat hessian_gamma(n_params_gamma, n_params_gamma);
    new_gamma.set_size(n_params_gamma);

    compute_mu(gamma_design_mat, prev_gamma, mu_gamma);
    compute_gradient(gamma_design_mat, w_t, Y_unval_all, mu_gamma, n, gradient_gamma);
    compute_hessian(gamma_design_mat, w_t, mu_gamma, n, hessian_gamma);

    gamma_success = newton_raphson_step(prev_gamma, gradient_gamma, hessian_gamma, new_gamma);

    arma::colvec gamma_diff = arma::abs(new_gamma - prev_gamma);
    gamma_conv = arma::all(gamma_diff < TOL);
  } else {
    new_gamma = prev_gamma;
  }

  // ===================================================================
  // Update p
  // ===================================================================
  arma::mat new_p(prev_p.n_rows, prev_p.n_cols);
  update_p(p_val_num, u_t, m, Nu, new_p);

  arma::mat p_diff = arma::abs(new_p - prev_p);
  bool p_conv = arma::all(arma::vectorise(p_diff) < TOL);

  // ===================================================================
  // Return results
  // ===================================================================
  return Rcpp::List::create(
    Rcpp::Named("theta") = new_theta,
    Rcpp::Named("gamma") = new_gamma,
    Rcpp::Named("p") = new_p,
    Rcpp::Named("theta_conv") = theta_conv,
    Rcpp::Named("gamma_conv") = gamma_conv,
    Rcpp::Named("p_conv") = p_conv,
    Rcpp::Named("theta_success") = theta_success,
    Rcpp::Named("gamma_success") = gamma_success
  );
}
