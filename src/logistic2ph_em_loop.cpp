#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Forward declarations of helper functions
inline void compute_mu(         const arma::mat&    design_mat,
                                const arma::colvec& beta,
                                arma::colvec&       mu_out);

inline void compute_gradient(   const arma::mat&    design_mat,
                                const arma::colvec& w_t,
                                const arma::colvec& y,
                                const arma::colvec& mu,
                                arma::colvec&       gradient_out);

inline void compute_hessian(    const arma::mat&    design_mat,
                                const arma::colvec& w_t,
                                const arma::colvec& mu,
                                arma::mat&          hessian_out);

inline bool newton_raphson_step(const arma::colvec& beta_old,
                                const arma::colvec& gradient,
                                const arma::mat&    hessian,
                                arma::colvec&       beta_new);

// Complete EM algorithm in C++
// [[Rcpp::export]]
Rcpp::List logistic2ph_em_loop(
    const arma::mat&    theta_design_mat,
    const arma::mat&    gamma_design_mat,
    const arma::mat&    comp_dat_all,
    const arma::mat&    comp_dat_unval,
    const arma::mat&    Bspline,
    const arma::colvec& Y_all,
    const arma::colvec& Y_unval_all,
    const arma::colvec& Y_unval_vec,
    const arma::colvec& Y_vec,
    const arma::mat&    p_val_num,
    arma::colvec        theta_init,
    arma::colvec        gamma_init,
    arma::mat           p_init,
    int                 n,
    int                 N,
    int                 m,
    bool                errorsY,
    double              TOL,
    int                 MAX_ITER,
    bool                verbose
)
{
  const int    Nu             = N - n;
  const int    n_spline       = Bspline.n_cols;
  const int    mult           = errorsY ? 2 : 1;
  const int    total_rows     = mult * m * Nu;
  const int    n_params_theta = theta_design_mat.n_cols;
  const int    n_params_gamma = gamma_design_mat.n_cols;
  const int    n_total        = theta_design_mat.n_rows;

  // ===================================================================
  // Pre-allocate ALL workspace memory
  // ===================================================================
  // E-step workspace
  arma::mat    psi_t(total_rows, n_spline);
  arma::colvec psi_denom(Nu);
  arma::colvec w_t_unval(total_rows);
  arma::mat    u_t(m * Nu, n_spline);
  arma::colvec w_t_lengthened(n_total);

  // M-step workspace
  arma::colvec mu_theta(n_total);
  arma::colvec gradient_theta(n_params_theta);
  arma::mat    hessian_theta(n_params_theta, n_params_theta);

  arma::colvec mu_gamma(n_total);
  arma::colvec gradient_gamma(n_params_gamma);
  arma::mat    hessian_gamma(n_params_gamma, n_params_gamma);

  arma::mat    new_p_num(m, n_spline);

  // Current parameter values
  arma::colvec theta      = theta_init;
  arma::colvec gamma      = gamma_init;
  arma::mat    p          = p_init;

  // Convergence tracking
  bool         converged  = false;
  int          iterations = 0;

  // ===================================================================
  // Main EM Loop
  // ===================================================================
  for (int it = 0; it < MAX_ITER; ++it)
  {
    iterations = it + 1;

    if (verbose && (it % 10 == 0 || it < 5))
    {
      Rcpp::Rcout << "Iteration: " << (it + 1) << std::endl;
    }

    // =================================================================
    // E-STEP
    // =================================================================
    // Use direct pointer access to first column as workspace
    double* py_ptr = psi_t.colptr(0);

    // Compute p(Y | X) directly in workspace
    arma::colvec mu_theta_estep = theta_design_mat.rows(n, theta_design_mat.n_rows - 1) * theta;
    for (int i = 0; i < total_rows; ++i)
    {
      py_ptr[i] = 1.0 / (1.0 + std::exp(-mu_theta_estep[i]));
      if (Y_vec[i] == 0) py_ptr[i] = 1.0 - py_ptr[i];
    }

    // Multiply by p(Y*|X*,Y,X) if needed (in-place)
    if (errorsY)
    {
      arma::colvec mu_gamma_estep = gamma_design_mat.rows(n, gamma_design_mat.n_rows - 1) * gamma;
      for (int i = 0; i < total_rows; ++i)
      {
        double p_val = 1.0 / (1.0 + std::exp(-mu_gamma_estep[i]));
        py_ptr[i] *= (Y_unval_vec[i] == 0) ? (1.0 - p_val) : p_val;
      }
    }

    // Compute psi_num directly in psi_t (no intermediate matrix)
    for (int y_idx = 0; y_idx < mult; ++y_idx)
    {
      int y_offset = y_idx * m * Nu;

      for (int k = 0; k < m; ++k)
      {
        int k_offset = y_offset + k * Nu;

        for (int i = 0; i < Nu; ++i)
        {
          int row_idx = k_offset + i;
          double py_val = py_ptr[row_idx];

          for (int j = 0; j < n_spline; ++j)
          {
            psi_t(row_idx, j) = py_val * p(k, j) * Bspline(row_idx, j);
          }
        }
      }
    }

    // Compute denominator
    psi_denom.zeros();
    for (int y_idx = 0; y_idx < mult; ++y_idx)
    {
      for (int k = 0; k < m; ++k)
      {
        int k_offset = y_idx * (m * Nu) + k * Nu;

        for (int i = 0; i < Nu; ++i)
        {
          int row_idx = k_offset + i;
          double row_sum = 0.0;
          for (int j = 0; j < n_spline; ++j)
          {
            row_sum += psi_t(row_idx, j);
          }
          psi_denom[i] += row_sum;
        }
      }
    }

    // Avoid division by zero
    for (int i = 0; i < Nu; ++i)
    {
      if (psi_denom[i] == 0.0) psi_denom[i] = 1.0;
    }

    // Normalize psi_t IN-PLACE
    for (int y_idx = 0; y_idx < mult; ++y_idx)
    {
      for (int k = 0; k < m; ++k)
      {
        int k_offset = y_idx * (m * Nu) + k * Nu;

        for (int i = 0; i < Nu; ++i)
        {
          int row_idx = k_offset + i;
          double denom = psi_denom[i];

          for (int j = 0; j < n_spline; ++j)
          {
            psi_t(row_idx, j) /= denom;
          }
        }
      }
    }

    // Compute w_t (sum across columns)
    for (int i = 0; i < total_rows; ++i)
    {
      double sum = 0.0;
      for (int j = 0; j < n_spline; ++j)
      {
        sum += psi_t(i, j);
      }
      w_t_unval[i] = sum;
    }

    // Compute u_t (collapse Y dimension)
    if (errorsY)
    {
      for (int k = 0; k < m; ++k)
      {
        for (int i = 0; i < Nu; ++i)
        {
          int out_row = k * Nu + i;
          int row_y0 = k * Nu + i;
          int row_y1 = (m * Nu) + k * Nu + i;

          for (int j = 0; j < n_spline; ++j)
          {
            u_t(out_row, j) = psi_t(row_y0, j) + psi_t(row_y1, j);
          }
        }
      }
    } else
    {
      u_t = psi_t.rows(0, m * Nu - 1);
    }

    // Lengthen w_t (prepend n ones)
    w_t_lengthened.head(n).ones();
    w_t_lengthened.tail(total_rows) = w_t_unval;

    // =================================================================
    // M-STEP
    // =================================================================
    arma::colvec theta_old = theta;
    arma::colvec gamma_old = gamma;
    arma::mat    p_old     = p;

    // Update theta
    compute_mu(theta_design_mat, theta, mu_theta);
    compute_gradient(theta_design_mat, w_t_lengthened, Y_all, mu_theta, gradient_theta);
    compute_hessian(theta_design_mat, w_t_lengthened, mu_theta, hessian_theta);

    bool theta_success = newton_raphson_step(theta, gradient_theta, hessian_theta, theta);

    if (!theta_success)
    {
      if (verbose) Rcpp::Rcout << "Warning: Hessian singular for theta at iteration " << (it + 1) << std::endl;
      return Rcpp::List::create(
        Rcpp::Named("converged") = false,
        Rcpp::Named("iterations") = iterations,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("gamma") = gamma,
        Rcpp::Named("p") = p,
        Rcpp::Named("theta_success") = false,
        Rcpp::Named("message") = "Singular Hessian for theta - fallback to R glm() needed"
      );
    }

    // Update gamma (if errorsY)
    if (errorsY)
    {
      compute_mu(gamma_design_mat, gamma, mu_gamma);
      compute_gradient(gamma_design_mat, w_t_lengthened, Y_unval_all, mu_gamma, gradient_gamma);
      compute_hessian(gamma_design_mat, w_t_lengthened, mu_gamma, hessian_gamma);

      bool gamma_success = newton_raphson_step(gamma, gradient_gamma, hessian_gamma, gamma);

      if (!gamma_success)
      {
        if (verbose) Rcpp::Rcout << "Warning: Hessian singular for gamma at iteration " << (it + 1) << std::endl;
        return Rcpp::List::create(
          Rcpp::Named("converged") = false,
          Rcpp::Named("iterations") = iterations,
          Rcpp::Named("theta") = theta,
          Rcpp::Named("gamma") = gamma,
          Rcpp::Named("p") = p,
          Rcpp::Named("gamma_success") = false,
          Rcpp::Named("message") = "Singular Hessian for gamma - fallback to R glm() needed"
        );
      }
    }

    // Update p
    new_p_num = p_val_num;
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

    arma::rowvec col_sums = arma::sum(new_p_num, 0);
    for (int k = 0; k < m; ++k)
    {
      for (int j = 0; j < n_spline; ++j)
      {
        p(k, j) = new_p_num(k, j) / col_sums[j];
      }
    }

    // =================================================================
    // Check convergence
    // =================================================================
    bool theta_conv = arma::all(arma::abs(theta - theta_old) < TOL);
    bool gamma_conv = !errorsY || arma::all(arma::abs(gamma - gamma_old) < TOL);
    bool p_conv = arma::all(arma::vectorise(arma::abs(p - p_old)) < TOL);

    if (theta_conv && gamma_conv && p_conv)
    {
      converged = true;
      if (verbose) Rcpp::Rcout << "Converged at iteration " << (it + 1) << std::endl;
      break;
    }
  }

  // ===================================================================
  // Return results
  // ===================================================================
  return Rcpp::List::create(
    Rcpp::Named("converged")     = converged,
    Rcpp::Named("iterations")    = iterations,
    Rcpp::Named("theta")         = theta,
    Rcpp::Named("gamma")         = gamma,
    Rcpp::Named("p")             = p,
    Rcpp::Named("theta_success") = true,
    Rcpp::Named("gamma_success") = true,
    Rcpp::Named("message")       = converged ? "Converged" : "MAX_ITER reached"
  );
}

// =====================================================================
// Helper function implementations
// =====================================================================

inline void compute_mu(
    const arma::mat&    design_mat,
    const arma::colvec& beta,
    arma::colvec&       mu_out)
{
  arma::colvec xb = design_mat * beta;
  xb              = -xb;
  mu_out          = arma::exp(xb) / (1.0 + arma::exp(xb));
}

inline void compute_gradient(
    const arma::mat&    design_mat,
    const arma::colvec& w_t,
    const arma::colvec& y,
    const arma::colvec& mu,
    arma::colvec&       gradient_out)
{
  arma::colvec sumsVector = y - 1.0 + mu;
  arma::mat temp          = design_mat;
  temp.each_col()        %= sumsVector;
  temp.each_col()        %= w_t;
  gradient_out            = arma::sum(temp, 0).t();
}

inline void compute_hessian(
    const arma::mat&    design_mat,
    const arma::colvec& w_t,
    const arma::colvec& mu,
    arma::mat&          hessian_out)
{
  arma::colvec mus            = mu % (mu - 1.0) % w_t;
  arma::mat weighted_design   = design_mat;
  weighted_design.each_col() %= mus;
  hessian_out                 = design_mat.t() * weighted_design;
}

inline bool newton_raphson_step(
    const arma::colvec& beta_old,
    const arma::colvec& gradient,
    const arma::mat&    hessian,
    arma::colvec&       beta_new)
{
  arma::colvec delta;
  bool         success = arma::solve(delta, hessian, -gradient, arma::solve_opts::likely_sympd);

  if (success)
  {
    beta_new = beta_old + delta;
    return true;
  } else
  {
    beta_new = beta_old;
    return false;
  }
}
