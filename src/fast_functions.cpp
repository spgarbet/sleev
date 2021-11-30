// This definition allows us to do some big matrix multiplication (calculateHessian)
// Set in Makevars
// #define ARMA_64BIT_WORD 1

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdlib.h>


using namespace Rcpp;
using namespace std;



arma::mat matTimesVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  // Never an issue in this project, but if you're copy+pasting this code into yours,
  // you may want to uncomment the error checking
  // if (v.n_elem < mat.n_rows)
  // {
  //   arma::vec oldV = v;
  //   for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
  //   {
  //     v = join_vert(v, oldV);
  //   }
  // }

  // Multiply each col by the vector
  mat.each_col() %= v;

  return mat;
}

arma::mat matDivideVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  if (v.n_elem < mat.n_rows)
  {
    arma::vec oldV = v;
    for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
    {
      v = join_vert(v, oldV);
    }
  }

  // Divide each col by the vector
  mat.each_col() /= v;

  return mat;
}


// TRANSLATING PACKAGE FUNCTIONS TO CPP FOR SPEED BOOST

arma::vec lengthenWT(
  const arma::vec& w_t_original,
  const int& n,
  const bool& modifyW_T = true)
{
  // If we don't modify w_t_original, return it immediately
  if (!modifyW_T)
    return w_t_original;

  // If w_t_original is empty, return the ones vector of length n
  if (w_t_original.is_empty())
      return arma::ones<arma::vec>(n);

  // Put n 1's in front of w_t_original
  arma::vec w_t(n + w_t_original.n_elem);   // Initialize vector
  w_t = join_vert(arma::ones<arma::vec>(n), w_t_original);  // Join vector of 1's with w_t_original


  return w_t;
}

arma::vec calculateMu(
  const arma::mat& design_mat,
  const arma::mat& prev)
{
  arma::mat mu = (design_mat * prev) * -1;
  arma::vec mu1 = exp(mu).as_col();
  return mu1 / (1 + mu1);
}

arma::vec calculateGradient(
  arma::vec& w_t,
  const int& n,
  const arma::mat& design_mat,
  const arma::vec& Y_col,
  const arma::vec& muVector,
  const bool& modifyW_T = false)
{
  // Put n 1's in front of w_t
  w_t = lengthenWT(w_t, n, modifyW_T);

  // Calculate gradient
  arma::vec sumsVector = Y_col - 1 + muVector;

  // Convert to NumericMatrix for Rcpp sugar's colSums
  arma::mat temp = matTimesVec(design_mat, sumsVector);
  temp = matTimesVec(temp, w_t);

  // Sum's default behavior on matrices is like colSums
  arma::rowvec gradient = sum(temp);


  return reshape(gradient, gradient.n_elem, 1);

}


arma::mat calculateHessian(
  const arma::mat& design_mat,
  arma::vec& w_t,
  const arma::vec& muVector,
  const int & n,
  arma::vec& mus,
  const bool& modifyW_T = false  )
{
  w_t = lengthenWT(w_t, n, modifyW_T);

  // post_multiply = c(w_t * muVector * (muVector - 1)) * gamma_design_mat
  mus = muVector % (muVector - 1);
  mus = mus % w_t;

  return design_mat.t() * matTimesVec(design_mat, mus);

}

arma::vec pYstarCalc(
  const arma::mat& gamma_design_mat,
  const int& n,
  const int& excludeRows,
  const arma::mat& prev_gamma,
  const arma::mat& comp_dat_all,
  const int& Y_unval_index,
  arma::vec& pYstar,
  arma::vec& mu_gamma  )
{
  // pYstar and mu_gamma are pre-allocated to save time with memory management

  // same as gamma_design_mat[-c(1:n),]
  // get the elements of gamma_design_mat excluding the first excludeRows rows
  arma::mat filtered_gamma_design_mat = gamma_design_mat.rows(excludeRows, gamma_design_mat.n_rows-1);

  mu_gamma = filtered_gamma_design_mat * prev_gamma;
  pYstar = 1 / (1 + exp(mu_gamma * -1));

  arma::vec checkVector = comp_dat_all.col(Y_unval_index).rows(n, comp_dat_all.n_rows-1);
  for (unsigned int i = 0; i < pYstar.size(); ++i)
  {
    if (checkVector(i) == 0)
    {
      pYstar(i) = 1 - pYstar(i);
    }
  }

  return pYstar;
}

arma::mat pXCalc(
  const int& n,
  const arma::mat& comp_dat_all_cropped,
  const bool& errorsX,
  const bool& errorsY,
  arma::mat& pX,
  arma::mat& prevRows)
{
  // pX and prevRows are pre-allocated
  // arma::mat prevRows = prev_p.rows(indices);

  if (errorsX and errorsY)
  {
    // need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    // multiply by the B-spline terms

    // these indices need to be repeated 2x
    arma::mat joinedPrevP = join_vert(prevRows, prevRows);

    // element-wise multiplication
    pX = joinedPrevP % comp_dat_all_cropped;
  }
  else if (errorsX)
  {
    // need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    // multiply by the B-spline terms
    pX = prevRows % comp_dat_all_cropped;
  }

  return pX;
}
