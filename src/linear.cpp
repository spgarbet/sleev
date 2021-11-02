#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "utility.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;


double WaldLinearGeneralSplineProfile (MatrixXd& pB, RowVectorXd& p_col_sum,
	VectorXd& q_row_sum, MatrixXd& p, MatrixXd& p0, MatrixXd& P_theta, MatrixXd& q, VectorXd& resi_n, MatrixXd& logp,
	const VectorXd& theta, const VectorXd& Y, const MatrixXd& X, const MatrixXd& Bspline_uni,
	const MatrixXd& ZW, const MatrixXd& X_uni, const VectorXi& X_uni_ind,
	const VectorXi& Bspline_uni_ind, const MatrixXd& p_static, const double& sigma_sq, const int& n, const int& n2, const int& m,
	const int& s, const int& n_minus_n2, const int& X_nc, const int& ZW_nc, const int& MAX_ITER, const double& TOL)
{
	/**** temporary variables **********************************************************************************************************************/
	double tol;
	int iter;
	/* test code */
	// time_t t1, t2;
	/* test code end */
	/**** temporary variables **********************************************************************************************************************/

	/**** update P_theta ***************************************************************************************************************************/
	P_theta.col(0) = Y.tail(n_minus_n2);
	P_theta.col(0).noalias() -= ZW.bottomRows(n_minus_n2)*theta.tail(ZW_nc);
	// P_theta.col(0) = Y.tail(n_minus_n2) - ZW.bottomRows(n_minus_n2) * theta.tail(ZW_nc);

	const VectorXd pThetaColZero = P_theta.col(0);
	for (int k = 0; k < m; ++k)
	{
		P_theta.col(k).noalias() = pThetaColZero - VectorXd::Constant(n_minus_n2, X_uni(k, 0) * theta(0));
	}
	P_theta = P_theta.array().square();
	P_theta /= -2.*sigma_sq;
	P_theta = P_theta.array().exp();
	/**** update P_theta ***************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	p_col_sum = p_static.colwise().sum();
	for (int j = 0; j < s; ++j)
	{
		p.col(j) = p_static.col(j) / p_col_sum(j);
	}
	p0 = p;
	MatrixXd Bspline_Mat(n_minus_n2, m);
 	MatrixXd pthetaOverQ(P_theta.rows(), P_theta.cols());

	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		/* test code */
		// auto time = tic();
		// time(&t1);
		/* test code end */

		/**** E-step *******************************************************************************************************************************/

		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
		/**** update pB ****************************************************************************************************************************/

		/**** update q, q_row_sum ******************************************************************************************************************/
		// Construct Bspline matrix - ~2x faster (0.01)

		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i + n2));
		}
		q = P_theta.array() * Bspline_Mat.array();

		// for (int i = 0; i < n_minus_n2; ++i)
		// {
		// 	for (int k = 0; k < m; ++k)
		// 	{
		// 		// Rcout << i << ", " << k << endl;
		// 		q(i,k) = P_theta(i,k) * pB(Bspline_uni_ind(i + n2), k);
		// 	}
		// }

		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
		/**** update q, q_row_sum ******************************************************************************************************************/

		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/

		/**** update p *****************************************************************************************************************************/
		p.setZero();
		// Each row of P_theta is divided by q_row_sum[row]
		pthetaOverQ = P_theta.array().colwise() / q_row_sum.array();

		for (int i = 0; i < n_minus_n2; ++i)
		{
			// for (int k = 0; k < m; ++k)
			// {
			// 	// for (int j=0; j<s; ++j)
			// 	// {
			// 		// p(k,j) += Bspline_uni(Bspline_uni_ind(i + n2), j) * P_theta(i,k) / q_row_sum(i);
			// 		 // p(k,j) += Bspline_uni(Bspline_uni_ind(i + n2), j) * pthetaOverQ(i*m + k);
			// 	// }
			// 	p.row(k) += Bspline_Mat.row(i) * pthetaOverQ(i,k);
			// }

			// Col vector * row vector outer product
			// Creates a matrix of size p made of every element of Bspline_mat * every element of pthetaOverQ
			// ~10 sec speedup
			p += pthetaOverQ.row(i).transpose() * Bspline_uni.row(Bspline_uni_ind(i + n2));
		}
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; ++j)
		{
			p.col(j) /= p_col_sum(j);
		}
		/**** update p *****************************************************************************************************************************/

		/**** M-step *******************************************************************************************************************************/

		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (p-p0).array().abs().sum();
		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

		/**** update parameters ********************************************************************************************************************/
		p0 = p;
		/**** update parameters ********************************************************************************************************************/

		/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
		/**** check convergence ********************************************************************************************************************/

		/* test code */
		// time(&t2);
		// Rcpp::Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// Rcout << iter << '\t' << chrono::duration<double>(tic() - time).count() << '\t' << tol << endl;
		/* test code end */
	}

	if (iter == MAX_ITER)
	{
		return -999.;
	}
	else
	{
		/**** calculate the likelihood *************************************************************************************************************/
		double tmp, loglik;

		logp = p.array().log();
		for (int k=0; k<m; ++k)
		{
			for (int j=0; j<s; ++j)
			{
				if (p(k,j) <= 0.)
				{
					logp(k,j) = 0.;
				}
			}
		}
		pB = Bspline_uni*logp.transpose();

		loglik = 0.;
		for (int i=0; i<n2; ++i)
		{
			loglik += pB(Bspline_uni_ind(i),X_uni_ind(i));
		}

		pB = Bspline_uni*p.transpose();
		Bspline_Mat.setZero();
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i + n2));
		}
		q = P_theta.array() * Bspline_Mat.array();


		// for (int i=0; i<n_minus_n2; ++i)
		// {
		// 	for (int k=0; k<m; ++k)
		// 	{
		// 		q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
		// 	}
		// }
		q_row_sum = q.rowwise().sum();

		loglik += q_row_sum.array().log().sum();
		loglik += -log(2.*M_PI*sigma_sq)*n/2.;

		resi_n = Y.head(n2);
		resi_n.noalias() -= X*theta.head(X_nc);
		resi_n.noalias() -= ZW.topRows(n2)*theta.tail(ZW_nc);

		tmp = resi_n.squaredNorm();
		tmp /= 2.*sigma_sq;
		loglik -= tmp;
		/**** calculate the likelihood *************************************************************************************************************/

		return loglik;
	}
} // WaldLinearGeneralSplineProfile



// [[Rcpp::export]]
List TwoPhase_GeneralSpline (
	const MapVecd& Y,
	const MapMatd& X,
	const MapMatd& ZW,
	const MapMatd& Bspline,
	const double& hn,
	const int& MAX_ITER,
	const double& TOL,
	const int& noSE)
{
	// auto start = tic();
	/*#############################################################################################################################################*/
	/**** pass arguments from R to cpp *************************************************************************************************************/
	// const MapVecd Y(as<MapVecd>(Y_R));
	// const MapMatd X(as<MapMatd>(X_R));
	// const MapMatd ZW(as<MapMatd>(ZW_R));
	// const MapMatd Bspline(as<MapMatd>(Bspline_R));
	// const double hn = NumericVector(hn_R)[0];
	// const int MAX_ITER = IntegerVector(MAX_ITER_R)[0];
	// const double TOL = NumericVector(TOL_R)[0];
	// const int noSE = IntegerVector(noSE_R)[0];
	/**** pass arguments from R to cpp *************************************************************************************************************/
	/*#############################################################################################################################################*/


	// auto time = tic();
	/*#############################################################################################################################################*/
	/**** some useful constants ********************************************************************************************************************/
	const int n = Y.size();  // number of subjects in the first phase
	const int n2 = X.rows(); // number of subjects in the second phase
	const int n_minus_n2 = n - n2; // number of subjects not selected in the second phase
	const int ZW_nc = ZW.cols(); // number of inexpensive covariates
	const int X_nc = X.cols(); // number of expensive covariates X
	const int ncov = X_nc + ZW_nc; // number of all covariates
	const int s = Bspline.cols(); // number of B-spline functions
	/**** some useful constants ********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	// m: number of distinct rows of X
	// X_uni_ind(n2): row index of rows of X in X_uni
	// X_uni(m, X_nc): distinct rows of X
	// X_uni_n(m): count of appearances of each distinct row of X
	VectorXi X_index = indexx_Matrix_Row(X);
	int m = Num_Uni_Matrix_Row(X, X_index);
	VectorXi X_uni_ind(n2);
	MatrixXd X_uni(m, X_nc);
	VectorXi X_uni_n(m);
	Create_Uni_Matrix_Row(X, X_index, X_uni, X_uni_ind, X_uni_n);
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** summarize observed distinct rows of Bspline **********************************************************************************************/
	// m_B: number of distinct rows of Bspline
	// Bspline_uni_ind(n): row index of rows of Bspline in Bspline_uni
	// Bspline_uni(m_B, s): distinct rows of Bspline
	// Bspline_uni_n(m_B): count of appearances of each distinct row of Bspline
	VectorXi Bspline_index = indexx_Matrix_Row(Bspline);
	int m_B = Num_Uni_Matrix_Row(Bspline, Bspline_index);
	VectorXi Bspline_uni_ind(n);
	MatrixXd Bspline_uni(m_B, s);
	VectorXi Bspline_uni_n(m_B);
	Create_Uni_Matrix_Row(Bspline, Bspline_index, Bspline_uni, Bspline_uni_ind, Bspline_uni_n);
	/**** summarize observed distinct rows of Bspline **********************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	// t(Y)*Y
	const double LS_YtY_static = Y.squaredNorm();

	// t(X,ZW)*Y
	VectorXd LS_XtY_static(ncov);
	LS_XtY_static.head(X_nc) = X.transpose() * Y.head(n2);
	LS_XtY_static.tail(ZW_nc) = ZW.transpose() * Y;

	// t(X,ZW) * (X,ZW)
	MatrixXd LS_XtX_static(ncov, ncov);
	LS_XtX_static.topLeftCorner(X_nc,X_nc) = X.transpose() * X;
	LS_XtX_static.topRightCorner(X_nc,ZW_nc) = X.transpose() * ZW.topRows(n2);
	LS_XtX_static.bottomRightCorner(ZW_nc,ZW_nc) = ZW.transpose() * ZW;

	// p
	MatrixXd p_static(m, s);
	p_static.setZero();
	for (int i = 0; i < n2; ++i)
	{
		p_static.row(X_uni_ind(i)) += Bspline.row(i);
	}
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); 			// regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; 		// flag of none convergence in the estimation of regression coefficients
	bool flag_nonconvergence_cov; 	// flag of none convergence in the estimation of covariance matrix
	double sigma_sq; 				// residual variance
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd LS_XtY(ncov);
	MatrixXd LS_XtX(ncov, ncov);
	VectorXd theta0(ncov);
	double sigma_sq0;
	MatrixXd p(m, s);
	MatrixXd p0(m, s);
	RowVectorXd p_col_sum(s);
	MatrixXd q(n_minus_n2, m);
	VectorXd q_row_sum(n_minus_n2);
	RowVectorXd q_col_sum(m);
	MatrixXd pB(m_B, m);
	MatrixXd P_theta(n_minus_n2, m);
	MatrixXd Bspline_Mat(n_minus_n2, m);
	double tol;
	int iter, idx;
	// /* RT's test code */
	// time_t t1, t2;
	// /* RT's test code end */
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/

	// Rcout << "initializing " << chrono::duration<double> (tic() - time).count() << endl;

	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	theta.setZero();
	theta0.setZero();
	sigma_sq = sigma_sq0 = Var(Y);

	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; ++j)
	{
		p.col(j) = p_static.col(j) / p_col_sum(j);
	}
	p0 = p;

	flag_nonconvergence = false;
	flag_nonconvergence_cov = false;
	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		// /* RT's test code */
		// time(&t1);
		// /* RT's test code end */

		/**** E-step *******************************************************************************************************************************/

		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
		/**** update pB ****************************************************************************************************************************/

		/**** update P_theta ***********************************************************************************************************************/
		P_theta.col(0) = Y.tail(n_minus_n2) - ZW.bottomRows(n_minus_n2)*theta.tail(ZW_nc);

		const VectorXd pThetaColZero = P_theta.col(0);
		for (int k=0; k<m; ++k)
		{
			P_theta.col(k).noalias() = pThetaColZero - VectorXd::Constant(n_minus_n2, X_uni(k,0) * theta(0));
		}
		P_theta = P_theta.array().square();
		P_theta /= -2. * sigma_sq;
		P_theta = P_theta.array().exp();
		/**** update P_theta ***********************************************************************************************************************/
		// time = tic();
		/**** update q, q_row_sum ******************************************************************************************************************/

		// for (int i = 0; i < n_minus_n2; ++i)
		// {
		// 	for (int k = 0; k < m; ++k)
		// 	{
		// 		q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
		// 	}
		// 	// q.row(i).array() = P_theta.row(i).array() * pB.row(Bspline_uni_ind(i+n2)).array();
		// }
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i+n2));
		}
		q = P_theta.array() * Bspline_Mat.array();

		q_row_sum = q.rowwise().sum();
		for (int i = 0; i < n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
		// q.array().colwise() /= q_row_sum;
		q_col_sum = q.colwise().sum();
		/**** update q, q_row_sum ******************************************************************************************************************/
		// Rcout << "update q " << chrono::duration<double> (tic() - time).count() << endl;
		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/
		// time= tic();
		/**** update theta and sigma_sq ************************************************************************************************************/
		LS_XtX = LS_XtX_static;
		LS_XtY = LS_XtY_static;

		for (int k=0; k<m; ++k)
		{
			LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += q_col_sum(k) * X_uni.row(k).transpose() * X_uni.row(k);
		}

		idx = n2;
		for (int i=0; i<n_minus_n2; ++i, ++idx)
		{
			for (int k=0; k<m; ++k)
			{
				const VectorXd firstProduct = q(i,k) * X_uni.row(k).transpose();
				LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += firstProduct * ZW.row(idx);
				LS_XtY.head(X_nc).noalias() += firstProduct * Y(idx);
			}
		}

		theta = LS_XtX.selfadjointView<Eigen::Upper>().ldlt().solve(LS_XtY);
		sigma_sq = LS_XtY.transpose()*theta;
		sigma_sq = (LS_YtY_static-sigma_sq)/n;
		/**** update theta and sigma_sq ************************************************************************************************************/
		// Rcout << "update theta & sigma_sq " << chrono::duration<double> (tic() - time).count() << endl;
		// time = tic();
		/**** update p *****************************************************************************************************************************/
		p.setZero();
		const MatrixXd pthetaOverQ = P_theta.array().colwise() / q_row_sum.array();
		for (int i=0; i<n_minus_n2; ++i)
		{
			// for (int k=0; k<m; ++k)
			// {
			// 	for (int j=0; j<s; ++j)
			// 	{
			// 		// O(n^4) !! mean 182 sec
			// 		p(k,j) += Bspline_uni(Bspline_uni_ind(i + n2), j) * P_theta(i,k) / q_row_sum(i);
			// 	}
			// 	// p.row(k) += Bspline_uni.row(Bspline_uni_ind(i+n2)) * pthetaOverQ(i,k);
			// }
			// mean 177 sec
			p += pthetaOverQ.row(i).transpose() * Bspline_uni.row(Bspline_uni_ind(i + n2));
		}
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j = 0; j < s; ++j)
		{
			p.col(j) /= p_col_sum(j);
		}
		/**** update p *****************************************************************************************************************************/
		// Rcout << "update p " << chrono::duration<double> (tic() - time).count() << endl;
		/**** M-step *******************************************************************************************************************************/


		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (theta-theta0).array().abs().sum();
		tol += fabs(sigma_sq-sigma_sq0);
		tol += (p-p0).array().abs().sum();
		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

		/**** update parameters ********************************************************************************************************************/
		theta0 = theta;
		sigma_sq0 = sigma_sq;
		p0 = p;
		/**** update parameters ********************************************************************************************************************/

		/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
		/**** check convergence ********************************************************************************************************************/

		// /* RT's test code */
		// time(&t2);
		// Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// /* RT's test code end */
	}
	/**** EM algorithm *****************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** variance estimation **********************************************************************************************************************/
	if (iter == MAX_ITER)
	{
		flag_nonconvergence = true;
		flag_nonconvergence_cov = true;
		theta.setConstant(-999.);
		sigma_sq = -999.;
		cov_theta.setConstant(-999.);
	}
	else if (noSE)
	{
		flag_nonconvergence_cov = true;
		cov_theta.setConstant(-999.);
	}
	else
	{
		VectorXd resi_n(n2);
		VectorXd profile_vec(ncov+1);
		MatrixXd logp(m, s);
		MatrixXd profile_mat(ncov+1, ncov+1);
		MatrixXd inv_profile_mat(ncov+1, ncov+1);
		double loglik;

		profile_mat.setZero();
		profile_vec.setZero();

		loglik = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n, logp, theta, Y, X, Bspline_uni, ZW, X_uni,
			X_uni_ind, Bspline_uni_ind, p_static, sigma_sq, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
		if (loglik == -999.)
		{
			flag_nonconvergence_cov = true;
		}
		// for (int i=0; i<ncov+1; ++i)
		// {
		// 	for (int j=i; j<ncov+1; ++j)
		// 	{
		// 		profile_mat(i,j) = loglik;
		// 	}
		// }
		profile_mat.triangularView<Upper>().setConstant(loglik);

		for (int i=0; i<ncov; ++i)
		{
			theta0 = theta;
			theta0(i) += hn;
			sigma_sq0 = sigma_sq;

			// The only thing that changes between iterations is theta0, which is barely used in the function
			// TODO: how to not rely on theta0 and just do the changes within the function

			profile_vec(i) = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n,
				logp, theta0, Y, X, Bspline_uni, ZW, X_uni,	X_uni_ind, Bspline_uni_ind, p_static, sigma_sq0,
				n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
		}
		theta0 = theta;
		sigma_sq0 = sigma_sq + hn;
		profile_vec(ncov) = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n,
		 logp, theta0, Y, X, Bspline_uni, ZW, X_uni, X_uni_ind, Bspline_uni_ind, p_static, sigma_sq0, n, n2, m,
		 s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
		// profile_vec = WaldLinearGeneralSplineProfileVector(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n,
		// 		logp, theta, Y, X, Bspline_uni, ZW, X_uni,	X_uni_ind, Bspline_uni_ind, p_static, sigma_sq,
		// 		n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL, hn, ncov);
		// for (int i = 0; i < ncov + 1; ++i)
		// {
		// 	Rcout << profile_vec(i) << endl;
		// 	if(profile_vec(i) == -999.)
		// 	{
		// 		flag_nonconvergence_cov = true;
		// 		break;
		// 	}
		// }
		flag_nonconvergence_cov = (profile_vec.array() == -999.).any();


		for (int i=0; i<ncov; ++i)
		{
			theta0 = theta;
			theta0(i) += hn;
			sigma_sq0 = sigma_sq+hn;
			loglik = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n, logp, theta0, Y, X, Bspline_uni, ZW, X_uni,
				X_uni_ind, Bspline_uni_ind, p_static, sigma_sq0, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
			Rcout << loglik << endl;
			if (loglik == -999.)
			{
				flag_nonconvergence_cov = true;
			}
			profile_mat(i,ncov) += loglik;
			profile_mat(i,ncov) -= profile_vec(i)+profile_vec(ncov);
			for (int j = i; j < ncov; ++j)
			{
				theta0 = theta;
				theta0(i) += hn;
				theta0(j) += hn;
				sigma_sq0 = sigma_sq;
				loglik = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n, logp, theta0, Y, X, Bspline_uni, ZW, X_uni,
					X_uni_ind, Bspline_uni_ind, p_static, sigma_sq0, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
				// Rcout << loglik << endl;
				if (loglik == -999.)
				{
					flag_nonconvergence_cov = true;
				}
				profile_mat(i,j) += loglik;
				profile_mat(i,j) -= profile_vec(i)+profile_vec(j);
			}
		}
		theta0 = theta;
		sigma_sq0 = sigma_sq+2.*hn;
		loglik = WaldLinearGeneralSplineProfile(pB, p_col_sum, q_row_sum, p, p0, P_theta, q, resi_n, logp, theta0, Y, X, Bspline_uni, ZW, X_uni,
			X_uni_ind, Bspline_uni_ind, p_static, sigma_sq0, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
		if (loglik == -999.)
		{
			flag_nonconvergence_cov = true;
		}
		profile_mat(ncov,ncov) += loglik;
		profile_mat(ncov,ncov) -= 2.*profile_vec(ncov);

		if (flag_nonconvergence_cov == true)
		{
			cov_theta.setConstant(-999.);
		}
		else
		{
			// for (int i=0; i<ncov+1; ++i)
			// {
			// 	for (int j=i+1; j<ncov+1; ++j)
			// 	{
			// 		profile_mat(j,i) = profile_mat(i,j);
			// 	}
			// }

			// BEFORE						AFTER
			// 0 	1 	2  	3 	4 			0 	1 	2  	3 	4
			// 5 	6 	7 	8 	9 			1 	6 	7 	8 	9
			// 10 	11 	12 	13 	14 			2 	7 	12 	13 	14
			// 15 	16 	17 	18 	19 			3 	8 	13 	18 	19
			// 20 	21 	22 	23 	24 			4 	9 	14 	19 	24
			profile_mat = profile_mat.selfadjointView<Upper>();

			profile_mat /= hn * hn;
			profile_mat = -profile_mat;
			// inv_profile_mat = profile_mat.selfadjointView<Upper>().ldlt().solve(MatrixXd::Identity(ncov+1, ncov+1));
			inv_profile_mat = profile_mat.ldlt().solve(MatrixXd::Identity(ncov+1, ncov+1));
			cov_theta = inv_profile_mat.topLeftCorner(ncov,ncov);

		}
	}
	/**** variance estimation **********************************************************************************************************************/
	/*#############################################################################################################################################*/


	// Rcout << "TwoPhase_GeneralSpline: " << chrono::duration<double> (tic() - start).count() << endl;

	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("theta") = theta,
		Named("sigma_sq") = sigma_sq,
		Named("cov_theta") = cov_theta,
		Named("flag_nonconvergence") = flag_nonconvergence,
		Named("flag_nonconvergence_cov") = flag_nonconvergence_cov);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // TwoPhase_GeneralSpline

MatrixXd WaldLinearVarianceMLE0 (const MatrixXd& LS_XtX, const VectorXd& LS_XtY, const VectorXd& p,
	const VectorXd& Y, const MatrixXd& ZW, const VectorXd& theta, const MatrixXd& X_uni, const MatrixXd& q,
	const int& n_minus_n2, const int& X_nc, const int& ZW_nc, const int& m, const int& ncov, const int& n, const int& n2, const double& sigma_sq)
{
	const int dimQ = ncov+1+m;

	VectorXd l1i(dimQ);
	VectorXd l1ik(dimQ);
	MatrixXd Q(dimQ, dimQ);
	MatrixXd resi(n_minus_n2, m);
	MatrixXd inv_profile_mat(dimQ-1, dimQ-1);
	// auto time = tic();
	/**** calculate resi ***************************************************************************************************************************/
	resi.col(0) = Y.tail(n_minus_n2);
	resi.col(0).noalias() -= ZW.bottomRows(n_minus_n2)*theta.tail(ZW_nc);

	for (int k=1; k<m; ++k)
	{
		resi.col(k) = resi.col(0);
	}

	for (int k=0; k<m; ++k)
	{
		resi.col(k).noalias() -= (X_uni.row(k)*theta.head(X_nc)).replicate(n_minus_n2,1);
	}
	/**** calculate resi ***************************************************************************************************************************/
	// Rcout << "calc resi: " << chrono::duration<double> (tic() - time).count() << endl;


	/**** augment the upper diagonal of Q **********************************************************************************************************/
	// time = tic();
	// add l2
	Q.setZero();
	Q.topLeftCorner(ncov,ncov) = LS_XtX/sigma_sq;
	Q.block(0,ncov,ncov,1) = (LS_XtY-LS_XtX*theta)/(sigma_sq*sigma_sq);
	Q(ncov,ncov) = (n+0.)/(2*sigma_sq*sigma_sq);
	for (int k=0; k<m; ++k)
	{
		Q(ncov+1+k,ncov+1+k) = (n+0.)/p(k);
	}
	// Rcout << "add 12: " << chrono::duration<double> (tic() - time).count() << endl;
	// time = tic();
	// add l1i, l1ik
	// add l1i: 1816.85 sec
	const double half = -1. / (2 * sigma_sq);
	const double sq = pow(2 * sigma_sq * sigma_sq, -1);
	const VectorXd head = X_uni / sigma_sq;
	const VectorXd oneOverP = 1. / p.array();
	const MatrixXd squaredResi = resi.array() * resi.array() * sq;
	for (int i = 0; i < n_minus_n2; ++i)
	{
		l1i.setZero();
		const VectorXd segment = (ZW.row(i+n2).transpose()) / sigma_sq;
		// auto innerloop = tic();
		for (int k = 0; k < m; ++k)
		{
			// auto one = tic();
			l1ik.setZero();
			// l1ik.head(X_nc) = resi(i,k) * (X_uni.row(k).transpose()) / sigma_sq;
			l1ik.head(X_nc) = resi(i,k) * head.row(k).transpose();
			l1ik.segment(X_nc,ZW_nc) = resi(i,k) * segment;
			// Rcout << "1: " << chrono::duration<double> (tic() - one).count() << endl;
			// one = tic();
			l1ik(ncov) = half + squaredResi(i,k) ;
			l1ik(ncov+1+k) = oneOverP(k);
			l1i += q(i,k) * l1ik;

			Q -= q(i,k) * l1ik * l1ik.transpose();
			// Rcout << q(i,k) * l1ik * l1ik.transpose() << endl;
			// Q -= q * l * l.t() (n_minus_n2 * m times)
			// Q += q * l * q * l.t() (n_minus_n2 times)
			// Rcout << "2: " << chrono::duration<double> (tic() - one).count() << endl;

		}


		Q += l1i * l1i.transpose();
		// Rcout << "inner loop: " << chrono::duration<double> (tic() - innerloop).count() << endl;
	}
	// Rcout << "add l1i: " << chrono::duration<double> (tic() - time).count() << endl;
	// time = tic();
	for (int k = 0; k < m-1; ++k)
	{
		Q.block(0, ncov+1+k, ncov+1, 1) -= Q.topRightCorner(ncov+1, 1);
		for (int kk = k; kk < m-1; ++kk)
		{
			Q(ncov+1+k, ncov+1+kk) -= Q(ncov+1+k, dimQ-1) + Q(ncov+1+kk, dimQ-1) - Q(dimQ-1, dimQ-1);
		}
	}
	// Rcout << "loop 2: " << chrono::duration<double> (tic() - time).count() << endl;
	/**** augment the upper diagonal of Q **********************************************************************************************************/
	inv_profile_mat = Q.topLeftCorner(dimQ-1,dimQ-1).selfadjointView<Eigen::Upper>().ldlt().solve(MatrixXd::Identity(dimQ-1, dimQ-1));
	return inv_profile_mat.topLeftCorner(ncov,ncov);

} // WaldLinearVarianceMLE0

// [[Rcpp::export]]
List TwoPhase_MLE0 (const VectorXd& Y, const MatrixXd& X,  const MatrixXd& ZW, const int& MAX_ITER, const double& TOL, const int& noSE)
/**** when there is no Z ***************************************************************************************************************************/
{
	// Rcout << "yes we're working" <<endl;
	// auto start = tic();
	/*#############################################################################################################################################*/
	/**** pass arguments from R to cpp *************************************************************************************************************/
	// const MapVecd Y(as<MapVecd>(Y_R));
	// const MapMatd X(as<MapMatd>(X_R));
	// const MapMatd ZW(as<MapMatd>(ZW_R));
	// const int MAX_ITER = IntegerVector(MAX_ITER_R)[0];
	// const double TOL = NumericVector(TOL_R)[0];
	// const int noSE = IntegerVector(noSE_R)[0];
	/**** pass arguments from R to cpp *************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some useful constants ********************************************************************************************************************/
	const int n = Y.size();  // number of subjects in the first phase
	const int n2 = X.rows(); // number of subjects in the second phase
	const int n_minus_n2 = n-n2; // number of subjects not selected in the second phase
	const int X_nc = X.cols(); // number of expensive covariates X
	const int ZW_nc = ZW.cols(); // number of inexpensive covariates
	const int ncov = X_nc+ZW_nc; // number of all covariates
	/**** some useful constants ********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	// m: number of distinct rows of X
	// X_uni_ind(n2): row index of rows of X in X_uni
	// X_uni(m, X_nc): distinct rows of X
	// X_uni_n(m): count of appearances of each distinct row of X
	VectorXi X_index = indexx_Matrix_Row(X);
	int m = Num_Uni_Matrix_Row(X, X_index);
	VectorXi X_uni_ind(n2);
	MatrixXd X_uni(m, X_nc);
	VectorXi X_uni_n(m);
	// modifies X_uni, X_uni_ind, and X_uni_n
	Create_Uni_Matrix_Row(X, X_index, X_uni, X_uni_ind, X_uni_n);
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	// t(Y)*Y
	const double LS_YtY_static = Y.squaredNorm();

	// t(X,ZW)*Y
	VectorXd LS_XtY_static(ncov);
	LS_XtY_static.head(X_nc) = X.transpose()*Y.head(n2);
	LS_XtY_static.tail(ZW_nc) = ZW.transpose()*Y;

	// t(X,ZW) * (X,ZW)
	MatrixXd LS_XtX_static(ncov, ncov);
	LS_XtX_static.topLeftCorner(X_nc,X_nc) = X.transpose()*X;
	LS_XtX_static.topRightCorner(X_nc,ZW_nc) = X.transpose()*ZW.topRows(n2);
	LS_XtX_static.bottomRightCorner(ZW_nc,ZW_nc) = ZW.transpose()*ZW;
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); // regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; // flag of none convergence in the estimation of regression coefficients
	bool flag_nonconvergence_cov; // flag of none convergence in the estimation of covariance matrix
	double sigma_sq; // residual variance
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd LS_XtY(ncov);
	MatrixXd LS_XtX(ncov, ncov);
	VectorXd theta0(ncov);
	double sigma_sq0;
	VectorXd p(m);
	VectorXd p0(m);
	MatrixXd q(n_minus_n2, m);
	VectorXd q_row_sum(n_minus_n2);
	RowVectorXd q_col_sum(m);
	MatrixXd P_theta(n_minus_n2, m);
	double tol;
	int iter, idx;
	// /* RT's test code */
	// time_t t1, t2;
	// /* RT's test code end */
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	theta.setZero();
	theta0.setZero();
	sigma_sq = sigma_sq0 = Var(Y);
	for (int k=0; k<m; ++k)
	{
		p(k) = (X_uni_n(k) + 0.) / (n2 + 0.);
	}
	p0 = p;
	flag_nonconvergence = false;
	flag_nonconvergence_cov = false;
	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		// /* RT's test code */
		// time(&t1);
		// /* RT's test code end */

		/**** E-step *******************************************************************************************************************************/
		// auto time = tic();
		/**** update P_theta ***********************************************************************************************************************/
		P_theta.col(0) = Y.tail(n_minus_n2);
		P_theta.col(0).noalias() -= ZW.bottomRows(n_minus_n2)*theta.tail(ZW_nc);
		for (int k=1; k<m; ++k)
		{
			P_theta.col(k) = P_theta.col(0);
		}
		for (int k=0; k<m; ++k)
		{
			P_theta.col(k).noalias() -= VectorXd::Constant(n_minus_n2, (X_uni.row(k)*theta.head(X_nc))(0,0));
		}
		P_theta = P_theta.array().square();
		P_theta /= -2.*sigma_sq;
		P_theta = P_theta.array().exp();
		/**** update P_theta ***********************************************************************************************************************/
		// Rcout << "update ptheta: " << chrono::duration<double> (tic() - time).count() << endl;
		// time = tic();
		/**** update q, q_row_sum ******************************************************************************************************************/
		for (int i=0; i<n_minus_n2; ++i)
		{
			// for (int k=0; k<m; ++k)
			// {
			// 	q(i,k) = P_theta(i,k)*p(k);
			// }
			q.row(i) = P_theta.row(i) * p;
		}
		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
		q_col_sum = q.colwise().sum();
		/**** update q, q_row_sum ******************************************************************************************************************/
		// Rcout << "update q: " << chrono::duration<double> (tic() - time).count() << endl;
		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/
		// time = tic();
		/**** update theta and sigma_sq ************************************************************************************************************/
		LS_XtX = LS_XtX_static;
		LS_XtY = LS_XtY_static;

		for (int k=0; k<m; ++k)
		{
			LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += q_col_sum(k)*X_uni.row(k).transpose()*X_uni.row(k);
		}

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			for (int k=0; k<m; ++k)
			{
				LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += q(i,k)*X_uni.row(k).transpose()*ZW.row(idx);
				LS_XtY.head(X_nc).noalias() += q(i,k)*X_uni.row(k).transpose()*Y(idx);
			}
		}

		theta = LS_XtX.selfadjointView<Eigen::Upper>().ldlt().solve(LS_XtY);
		sigma_sq = LS_XtY.transpose()*theta;
		sigma_sq = (LS_YtY_static-sigma_sq)/n;
		/**** update theta and sigma_sq ************************************************************************************************************/
		// Rcout << "update theta & sigma_sq: " << chrono::duration<double> (tic() - time).count() << endl;
		// time = tic();
		/**** update p *****************************************************************************************************************************/
		for (int k=0; k<m; ++k)
		{
			p(k) = X_uni_n(k)+0.;
		}
		p += q_col_sum.transpose();
		p /= n+0.;
		/**** update p *****************************************************************************************************************************/
		// Rcout << "update p: " << chrono::duration<double> (tic() - time).count() << endl;
		/**** M-step *******************************************************************************************************************************/


		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (theta-theta0).array().abs().sum();
		tol += fabs(sigma_sq-sigma_sq0);
		tol += (p-p0).array().abs().sum();
		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

		/**** update parameters ********************************************************************************************************************/
		theta0 = theta;
		sigma_sq0 = sigma_sq;
		p0 = p;
		/**** update parameters ********************************************************************************************************************/

		/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
		/**** check convergence ********************************************************************************************************************/

		// /* RT's test code */
		// time(&t2);
		// Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// /* RT's test code end */
	}
	/**** EM algorithm *****************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** variance estimation **********************************************************************************************************************/
	if (iter == MAX_ITER)
	{
		flag_nonconvergence = true;
		flag_nonconvergence_cov = true;
		theta.setConstant(-999.);
		sigma_sq = -999.;
		cov_theta.setConstant(-999.);
	}
	else if (noSE)
	{
		flag_nonconvergence_cov = true;
		cov_theta.setConstant(-999.);
	}
	else
	{
		cov_theta = WaldLinearVarianceMLE0(LS_XtX, LS_XtY, p, Y, ZW, theta, X_uni, q, n_minus_n2, X_nc, ZW_nc, m, ncov, n, n2, sigma_sq);
	}
	/**** variance estimation **********************************************************************************************************************/
	/*#############################################################################################################################################*/

	// Rcout << "TwoPhase_MLE0: " << chrono::duration<double> (tic() - start).count() << endl;


	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("theta") = theta,
		Named("sigma_sq") = sigma_sq,
		Named("cov_theta") = cov_theta,
		Named("flag_nonconvergence") = flag_nonconvergence,
		Named("flag_nonconvergence_cov") = flag_nonconvergence_cov);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // TwoPhase_MLE0
