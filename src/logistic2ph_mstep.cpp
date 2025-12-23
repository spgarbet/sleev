#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function: compute mu = exp(-X*beta) / (1 + exp(-X*beta))
// Matches .calculateMu from R
inline void compute_mu(
    const arma::mat& design_mat,
    const arma::colvec& beta,
    arma::colvec& mu_out)
{
  arma::colvec xb = design_mat * beta;
  xb = -xb;  // Negate
  mu_out = arma::exp(xb) / (1.0 + arma::exp(xb));
}

// Compute gradient: sum_i w_i * (y_i - 1 + mu_i) * x_i
// Matches .calculateGradient from R
// NOTE: w_t should already have ones prepended (lengthened)
inline void compute_gradient(
    const arma::mat& design_mat,
    const arma::colvec& w_t,
    const arma::colvec& y,
    const arma::colvec& mu,
    arma::colvec& gradient_out)
{
  int n_total = design_mat.n_rows;
  int n_params = design_mat.n_cols;

  // sumsVector = Y - 1 + mu
  arma::colvec sumsVector = y - 1.0 + mu;

  // Multiply design_mat by sumsVector element-wise per row, then by w_t
  // This matches: matTimesVec(matTimesVec(design_mat, sumsVector), w_t)
  arma::mat temp = design_mat;
  temp.each_col() %= sumsVector;  // Multiply each column by sumsVector
  temp.each_col() %= w_t;          // Multiply each column by w_t

  // Sum columns (like colSums in R)
  gradient_out = arma::sum(temp, 0).t();
}

// Compute Hessian: design_mat^T * diag(w_t * mu * (mu - 1)) * design_mat
// Matches .calculateHessian from R
// NOTE: The R version uses mu * (mu - 1), which is NEGATIVE
inline void compute_hessian(
    const arma::mat& design_mat,
    const arma::colvec& w_t,
    const arma::colvec& mu,
    arma::mat& hessian_out)
{
  int n_params = design_mat.n_cols;

  // mus = mu * (mu - 1) * w_t
  arma::colvec mus = mu % (mu - 1.0) % w_t;

  // Hessian = design_mat^T * diag(mus) * design_mat
  // Efficiently: design_mat^T * (design_mat with each col multiplied by mus)
  arma::mat weighted_design = design_mat;
  weighted_design.each_col() %= mus;

  hessian_out = design_mat.t() * weighted_design;
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

  for (int k = 0; k < m; ++k)
  {
    for (int i = 0; i < Nu; ++i)
    {
      int row_idx = k * Nu + i;
      for (int j = 0; j < n_spline; ++j)
      {
        new_p_num(k, j) += u_t(row_idx, j);
      }
    }
  }

  // Normalize: p_new = new_p_num / colSums(new_p_num)
  arma::rowvec col_sums = arma::sum(new_p_num, 0);

  for (int k = 0; k < m; ++k)
  {
    for (int j = 0; j < n_spline; ++j)
    {
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

  // NOTE: w_t should already be lengthened with n ones prepended by caller

  // ===================================================================
  // Update theta
  // ===================================================================
  arma::colvec mu_theta(theta_design_mat.n_rows);
  arma::colvec gradient_theta(n_params_theta);
  arma::mat hessian_theta(n_params_theta, n_params_theta);
  arma::colvec new_theta(n_params_theta);

  compute_mu(theta_design_mat, prev_theta, mu_theta);
  compute_gradient(theta_design_mat, w_t, Y_all, mu_theta, gradient_theta);
  compute_hessian(theta_design_mat, w_t, mu_theta, hessian_theta);

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
    compute_gradient(gamma_design_mat, w_t, Y_unval_all, mu_gamma, gradient_gamma);
    compute_hessian(gamma_design_mat, w_t, mu_gamma, hessian_gamma);

    gamma_success = newton_raphson_step(prev_gamma, gradient_gamma, hessian_gamma, new_gamma);

    arma::colvec gamma_diff = arma::abs(new_gamma - prev_gamma);
    gamma_conv = arma::all(gamma_diff < TOL);
  } else
  {
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
