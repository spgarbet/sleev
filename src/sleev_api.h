#ifndef SLEEV_API_H
#define SLEEV_API_H

#include <RcppArmadillo.h>

arma::vec pYstarCalc(
    const arma::mat& gamma_design_mat,
    const int&       startRow,
    const arma::mat& prev_gamma,
    const arma::mat& comp_dat_all,
    const int&       Y_unval_index
);

arma::vec lengthenWT(
    const arma::vec& w_t_original,
    const int&       n,
    const bool&      modifyW_T);

arma::vec calculateMu(
    const arma::mat& design_mat,
    const arma::mat& prev
);

arma::vec calculateGradient(
    arma::vec&       w_t,
    const int&       n,
    const arma::mat& design_mat,
    const arma::vec& Y_col,
    const arma::vec& muVector,
    const bool&      modifyW_T
);

arma::mat calculateHessian(
    const arma::mat& design_mat,
    arma::vec&       w_t,
    const arma::vec& muVector,
    const int& n,    arma::vec& mus,
    const bool&      modifyW_T
);

#endif
