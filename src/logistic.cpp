#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "utility.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;




double WaldLogisticGeneralSplineProfile (MatrixXd& pB, RowVectorXd& p_col_sum, VectorXd& ZW_theta, VectorXd& X_uni_theta,
	MatrixXd& e_X_uni_theta, MatrixXd& logi_X_uni_theta, VectorXd& e_X_theta,
	VectorXd& q_row_sum, MatrixXd& p, MatrixXd& p0, MatrixXd& P_theta, MatrixXd& q, MatrixXd& logp,
	const VectorXd& theta, const VectorXd& Y, const MatrixXd& X, const MatrixXd& Bspline_uni,
	const MatrixXd& ZW, const MatrixXd& X_uni, const VectorXi& X_uni_ind,
	const VectorXi& Bspline_uni_ind, const MatrixXd& p_static, const int n, const int n2, const int m,
	const int s, const int n_minus_n2, const int X_nc, const int ZW_nc, const int MAX_ITER, const double TOL)
{
	/**** temporary variables **********************************************************************************************************************/
	double tol;
	int iter, idx;
	/* test code */
	// time_t t1, t2;
	/* test code end */
	/**** temporary variables **********************************************************************************************************************/
	auto time = tic();
	/**** update P_theta ***************************************************************************************************************************/
	ZW_theta = ZW*theta.tail(ZW_nc);
	X_uni_theta = X_uni*theta.head(X_nc);
	for (int i=0; i<n_minus_n2; ++i)
	{
		for (int k=0; k<m; ++k)
		{
			e_X_uni_theta(i,k) = ZW_theta(i+n2)+X_uni_theta(k);
		}
	}
	e_X_uni_theta = e_X_uni_theta.array().exp();
	logi_X_uni_theta = e_X_uni_theta.array()/(e_X_uni_theta.array()+1.);

	for (int i=0; i<n_minus_n2; ++i)
	{
		idx = i+n2;
		if (Y(idx) == 1.)
		{
			P_theta.row(i) = logi_X_uni_theta.row(i);
		}
		else if (Y(idx) == 0.)
		{
			P_theta.row(i) = 1. - logi_X_uni_theta.row(i).array();
		}
	}
	/**** update P_theta ***************************************************************************************************************************/
	Rcout << "update p_theta " << chrono::duration<double>(tic() - time).count() << endl;
	/**** parameter initialization *****************************************************************************************************************/
	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; ++j)
	{
		p.col(j) = p_static.col(j)/p_col_sum(j);
	}
	p0 = p;
	MatrixXd Bspline_Mat(n_minus_n2, pB.cols());
	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		/* test code */
		// time(&t1);
		/* test code end */

		/**** E-step *******************************************************************************************************************************/

		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni * p.transpose();
		/**** update pB ****************************************************************************************************************************/
		time = tic();
		/**** update q, q_row_sum ******************************************************************************************************************/
		for (int i = 0; i < n_minus_n2; ++i)
		{
			Bspline_Mat.row(i) = pB.row(Bspline_uni_ind(i + n2));
		}
		q = P_theta.array() * Bspline_Mat.array();
		// for (int i=0; i<n_minus_n2; ++i)
		// {
		// 	// for (int k=0; k<m; ++k)
		// 	// {
		// 	// 	q(i,k) = P_theta(i,k) * pB(Bspline_uni_ind(i+n2), k);
		// 	// }
		// }
		// Rcout << q.cols() << endl << m << endl << P_theta.cols() << endl << pB.cols() << endl;

		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; ++i)
		{
			q.row(i) /= q_row_sum(i);
		}
		/**** update q, q_row_sum ******************************************************************************************************************/
		Rcout << "update q " << chrono::duration<double>(tic() - time).count() << endl;
		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/
		time = tic();
		/**** update p *****************************************************************************************************************************/
		p.setZero();

		// Each row of P_theta is divided by q_row_sum[row]
		const MatrixXd pthetaOverQ = P_theta.array().colwise() / q_row_sum.array();

		for (int i=0; i<n_minus_n2; ++i)
		{
			// for (int k=0; k<m; ++k)
			// {
			// 	for (int j=0; j<s; ++j)
			// 	{
			// 		p(k,j) += Bspline_uni(Bspline_uni_ind(i+n2),j)*P_theta(i,k)/q_row_sum(i);
			// 	}
			// }
			p += pthetaOverQ.row(i).transpose() * Bspline_uni.row(Bspline_uni_ind(i + n2));

			Rcpp::checkUserInterrupt();

		}

		p = p.array() * p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; ++j)
		{
			p.col(j) /= p_col_sum(j);
		}
		/**** update p *****************************************************************************************************************************/
		Rcout << "update p " << chrono::duration<double>(tic() - time).count() << endl;
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
		/* test code end */
	}

	if (iter == MAX_ITER)
	{
		return -999.;
	}
	else
	{
		/**** calculate the likelihood *************************************************************************************************************/
		double loglik;
		time = tic();

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
		for (int i=0; i<n2; ++i) {
			loglik += pB(Bspline_uni_ind(i),X_uni_ind(i));
		}

		pB = Bspline_uni*p.transpose();;

		// for (int i=0; i<n_minus_n2; ++i)
		// {	//problem
			// for (int k=0; k<m; ++k)
			// {
			// 	q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
			// }
			// q.row(i) = P_theta.row(i) * pB.row(Bspline_uni_ind(i+n2));
			// Rcpp::checkUserInterrupt();
		// }
		// Reusing Bspline_Mat
		q = P_theta.array() * Bspline_Mat.array();
		q_row_sum = q.rowwise().sum();

		loglik += q_row_sum.array().log().sum();

		for (int i=0; i<n2; ++i)
		{
			e_X_theta(i) = ZW_theta(i)+X_uni_theta(X_uni_ind(i));
			if (Y(i) == 1.)
			{
				loglik += e_X_theta(i);
			}
		}
		e_X_theta = e_X_theta.array().exp();
        e_X_theta = e_X_theta.array()+1.;
        loglik -= e_X_theta.array().log().sum();
		/**** calculate the likelihood *************************************************************************************************************/
        Rcout << "caclulate log likelihood " << chrono::duration<double>(tic() - time).count() << endl;
		return loglik;
	}
} // WaldLogisticGeneralSplineProfile

// [[Rcpp::export]]
List TwoPhase_GeneralSpline_logistic (const VectorXd& Y, const MatrixXd& X, const MatrixXd& ZW, const MatrixXd& Bspline, const double& hn, const int& MAX_ITER, const double& TOL, const int& noSE)
{
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



	/*#############################################################################################################################################*/
	/**** some useful constants ********************************************************************************************************************/
	const int n = Y.size();  // number of subjects in the first phase
	const int n2 = X.rows(); // number of subjects in the second phase
	const int n_minus_n2 = n-n2; // number of subjects not selected in the second phase
	const int ZW_nc = ZW.cols(); // number of inexpensive covariates
	const int X_nc = X.cols(); // number of expensive covariates X
	const int ncov = X_nc+ZW_nc; // number of all covariates
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
	VectorXi Bspline_uni_ind(n);
	VectorXi Bspline_index = indexx_Matrix_Row(Bspline);
	int m_B = Num_Uni_Matrix_Row(Bspline, Bspline_index);
	MatrixXd Bspline_uni(m_B, s);
	VectorXi Bspline_uni_n(m_B);
	Create_Uni_Matrix_Row(Bspline, Bspline_index, Bspline_uni, Bspline_uni_ind, Bspline_uni_n);
	/**** summarize observed distinct rows of Bspline **********************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	// p
	MatrixXd p_static(m, s);
	p_static.setZero();
	for (int i=0; i<n2; ++i)
	{
		p_static.row(X_uni_ind(i)) += Bspline.row(i);
	}
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); // regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; // flag of none convergence in the estimation of regression coefficients
	bool flag_nonconvergence_cov; // flag of none convergence in the estimation of covariance matrix
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd LS_XtY(ncov);
	MatrixXd LS_XtX(ncov, ncov);
	VectorXd theta0(ncov);
	MatrixXd p(m, s);
	MatrixXd p0(m, s);
	RowVectorXd p_col_sum(s);
	MatrixXd q(n_minus_n2, m);
	VectorXd q_row_sum(n_minus_n2);
	MatrixXd pB(m_B, m);
	MatrixXd P_theta(n_minus_n2, m);
	double tol;
	int iter, idx;
	VectorXd ZW_theta(n);
	VectorXd X_uni_theta(m);
	MatrixXd e_X_uni_theta(n_minus_n2, m);
	MatrixXd logi_X_uni_theta(n_minus_n2, m);
	MatrixXd logi2_X_uni_theta(n_minus_n2, m);
	VectorXd e_X_theta(n2);
	VectorXd logi_X_theta(n2);
	VectorXd logi2_X_theta(n2);
	string error_message;
	// /* RT test block */
	// time_t t1, t2;
	// /* RT test block */
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	theta.setZero();
	theta0.setZero();

	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; ++j)
	{
		p.col(j) = p_static.col(j)/p_col_sum(j);
	}
	p0 = p;

	flag_nonconvergence = false;
	flag_nonconvergence_cov = false;
	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		// /* RT test block */
		//	time(&t1);
		//	/* RT test block */

		/**** E-step *******************************************************************************************************************************/

		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
		/**** update pB ****************************************************************************************************************************/

		/**** update P_theta ***********************************************************************************************************************/
		auto time = tic();
		ZW_theta = ZW*theta.tail(ZW_nc);
		X_uni_theta = X_uni*theta.head(X_nc);
		for (int i=0; i<n_minus_n2; i++)
		{
			for (int k=0; k<m; k++)
			{
				e_X_uni_theta(i,k) = ZW_theta(i+n2)+X_uni_theta(k);
			}
		}
		e_X_uni_theta = e_X_uni_theta.array().exp();
        logi_X_uni_theta = e_X_uni_theta.array()/(e_X_uni_theta.array()+1.);
        logi2_X_uni_theta = logi_X_uni_theta.array()/(e_X_uni_theta.array()+1.);

		for (int i=0; i<n_minus_n2; i++)
		{
			idx = i+n2;
			if (Y(idx) == 1.)
			{
				P_theta.row(i) = logi_X_uni_theta.row(i);
			}
			else if (Y(idx) == 0.)
			{
				P_theta.row(i) = 1.-logi_X_uni_theta.row(i).array();
			}
			else
			{
				std::ostringstream error_message;
				error_message << "Error: the " << i+n2 << "th value of Y is not 1 or 0!\n";
				stop(error_message.str());
			}
		}
		Rcout << "update p theta " << chrono::duration<double>(tic() - time).count() << endl;

		/**** update P_theta ***********************************************************************************************************************/

		/**** update q, q_row_sum ******************************************************************************************************************/
		time = tic();
		for (int i=0; i<n_minus_n2; i++)
		{
			for (int k=0; k<m; k++)
			{
				q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
			}
		}
		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; i++)
		{
			q.row(i) /= q_row_sum(i);
		}
		Rcout << "update q's " << chrono::duration<double>(tic() - time).count() << endl;
		/**** update q, q_row_sum ******************************************************************************************************************/

		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/

		/**** update theta *************************************************************************************************************************/
		time = tic();
		for (int i=0; i<n2; i++)
		{
			e_X_theta(i) = ZW_theta(i)+X_uni_theta(X_uni_ind(i));
		}
		e_X_theta = e_X_theta.array().exp();
        logi_X_theta = e_X_theta.array()/(e_X_theta.array()+1.);
        logi2_X_theta = logi_X_theta.array()/(e_X_theta.array()+1.);

		LS_XtX.setZero();
		LS_XtY.setZero();
		for (int i=0; i<n2; i++)
		{
			LS_XtY.head(X_nc).noalias() += X.row(i).transpose()*(Y(i)-logi_X_theta(i));
			LS_XtY.tail(ZW_nc).noalias() += ZW.row(i).transpose()*(Y(i)-logi_X_theta(i));

			LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += (X.row(i).transpose())*X.row(i)*logi2_X_theta(i);
			LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += (X.row(i).transpose())*ZW.row(i)*logi2_X_theta(i);
			LS_XtX.bottomRightCorner(ZW_nc,ZW_nc).noalias() += (ZW.row(i).transpose())*ZW.row(i)*logi2_X_theta(i);
		}
		for (int i=0; i<n_minus_n2; i++)
		{
			idx = i+n2;
			for (int k=0; k<m; k++)
			{
				LS_XtY.head(X_nc).noalias() += q(i,k)*X_uni.row(k).transpose()*(Y(idx)-logi_X_uni_theta(i,k));
				LS_XtY.tail(ZW_nc).noalias() += q(i,k)*ZW.row(idx).transpose()*(Y(idx)-logi_X_uni_theta(i,k));

				LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += q(i,k)*(X_uni.row(k).transpose())*X_uni.row(k)*logi2_X_uni_theta(i,k);
				LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += q(i,k)*(X_uni.row(k).transpose())*ZW.row(idx)*logi2_X_uni_theta(i,k);
				LS_XtX.bottomRightCorner(ZW_nc,ZW_nc).noalias() += q(i,k)*(ZW.row(idx).transpose())*ZW.row(idx)*logi2_X_uni_theta(i,k);
			}
		}

		theta = LS_XtX.selfadjointView<Eigen::Upper>().ldlt().solve(LS_XtY);
		theta += theta0;
		Rcout << "update theta " << chrono::duration<double>(tic() - time).count() << endl;
		/**** update theta *************************************************************************************************************************/

		/**** update p *****************************************************************************************************************************/
		time = tic();
		p.setZero();
		for (int i=0; i<n_minus_n2; i++)
		{
			for (int k=0; k<m; k++)
			{
				for (int j=0; j<s; j++)
				{
					p(k,j) += Bspline_uni(Bspline_uni_ind(i+n2),j)*P_theta(i,k)/q_row_sum(i);
				}
			}
		}
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; j++)
		{
			p.col(j) /= p_col_sum(j);
		}
		Rcout << "update p" << chrono::duration<double>(tic() - time).count() << endl;
		/**** update p *****************************************************************************************************************************/

		/**** M-step *******************************************************************************************************************************/


		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (theta-theta0).array().abs().sum();
		tol += (p-p0).array().abs().sum();
		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

		/**** update parameters ********************************************************************************************************************/
		theta0 = theta;
		p0 = p;
		/**** update parameters ********************************************************************************************************************/

		/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
		/**** check convergence ********************************************************************************************************************/

		// /* RT test block */
		//	time(&t2);
		//	Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// /* RT test block */
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
		cov_theta.setConstant(-999.);
	}
	else if (noSE)
	{
		flag_nonconvergence_cov = true;
		cov_theta.setConstant(-999.);
	}
	else
	{
		VectorXd profile_vec(ncov);
		MatrixXd logp(m, s);
		MatrixXd profile_mat(ncov, ncov);
		double loglik;

		profile_mat.setZero();
		profile_vec.setZero();

		loglik = WaldLogisticGeneralSplineProfile(pB, p_col_sum, ZW_theta, X_uni_theta, e_X_uni_theta, logi_X_uni_theta, e_X_theta, q_row_sum, p, p0, P_theta,
			q, logp, theta, Y, X, Bspline_uni, ZW, X_uni, X_uni_ind, Bspline_uni_ind, p_static, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
		if (loglik == -999.)
		{
			flag_nonconvergence_cov = true;
		}
		// for (int i=0; i<ncov; i++)
		// {
		// 	for (int j=i; j<ncov; j++)
		// 	{
		// 		profile_mat(i,j) = loglik;
		// 	}
		// }
		profile_mat.triangularView<Upper>().setConstant(loglik);

		for (int i=0; i<ncov; ++i)
		{
			theta0 = theta;
			theta0(i) += hn;
			profile_vec(i) = WaldLogisticGeneralSplineProfile(pB, p_col_sum, ZW_theta, X_uni_theta, e_X_uni_theta, logi_X_uni_theta, e_X_theta, q_row_sum, p, p0, P_theta,
				q, logp, theta0, Y, X, Bspline_uni, ZW, X_uni, X_uni_ind, Bspline_uni_ind, p_static, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
			if(profile_vec(i) == -999.)
			{
				flag_nonconvergence_cov = true;
			}
		}

		for (int i=0; i<ncov; ++i)
		{
			for (int j=i; j<ncov; ++j)
			{
				theta0 = theta;
				theta0(i) += hn;
				theta0(j) += hn;
				loglik = WaldLogisticGeneralSplineProfile(pB, p_col_sum, ZW_theta, X_uni_theta, e_X_uni_theta, logi_X_uni_theta, e_X_theta, q_row_sum, p, p0, P_theta,
					q, logp, theta0, Y, X, Bspline_uni, ZW, X_uni, X_uni_ind, Bspline_uni_ind, p_static, n, n2, m, s, n_minus_n2, X_nc, ZW_nc, MAX_ITER, TOL);
				if (loglik == -999.)
				{
					flag_nonconvergence_cov = true;
				}
				profile_mat(i,j) += loglik;
				profile_mat(i,j) -= profile_vec(i)+profile_vec(j);
			}
		}

		if (flag_nonconvergence_cov == true)
		{
			cov_theta.setConstant(-999.);
		}
		else
		{
			// for (int i=0; i<ncov; i++)
			// {
			// 	for (int j=i+1; j<ncov; j++)
			// 	{
			// 		profile_mat(j,i) = profile_mat(i,j);
			// 	}
			// }
			profile_mat = profile_mat.selfadjointView<Upper>();

			profile_mat /= hn * hn;
			profile_mat = -profile_mat;
			cov_theta = profile_mat.selfadjointView<Upper>().ldlt().solve(MatrixXd::Identity(ncov, ncov));
			// TODO: is this the correct value for Identity matrix? In linear, it's ncov+1
		}
	}
	/**** variance estimation **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("theta") = theta,
						Named("cov_theta") = cov_theta,
						Named("flag_nonconvergence") = flag_nonconvergence,
						Named("flag_nonconvergence_cov") = flag_nonconvergence_cov);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // TwoPhase_GeneralSpline_logistic

MatrixXd WaldLogisticVarianceMLE0 (const MatrixXd& LS_XtX, const VectorXd& p, const MatrixXd& logi_X_uni_theta,
	const VectorXd& Y, const MatrixXd& ZW, const MatrixXd& X_uni, const MatrixXd& q,
	const int& n_minus_n2, const int& X_nc, const int& ZW_nc, const int& m, const int& ncov, const int& n, const int& n2)
{
	const int dimQ = ncov+m;
	const VectorXd oneOverP = p.pow(-1);

	VectorXd l1i(dimQ);
	VectorXd l1ik(dimQ);
	MatrixXd Q(dimQ, dimQ);
	MatrixXd resi(n_minus_n2, m);
	MatrixXd inv_profile_mat(dimQ-1, dimQ-1);
	int idx;

	/**** augment the upper diagonal of Q **********************************************************************************************************/
	// add l2
	Q.setZero();
	Q.topLeftCorner(ncov,ncov) = LS_XtX;
	for (int k=0; k<m; ++k) {
		Q(ncov + k, ncov + k) = (double)n * oneOverP(k);
	}

	// add l1i, l1ik
	MatrixXd yMinusLogi(n_minus_n2, m);
	for (int i = 0, idx = n2; i < n_minus_n2; ++i, ++idx)
	{
		yMinusLogi.row(i).setConstant(Y(idx));
	  	yMinusLogi.row(i) -= logi_X_uni_theta.row(i);
	}

	for (int i = 0, idx = n2; i < n_minus_n2; ++i, ++idx)
	// for (int i = 0, idx = n2; i < 100; ++i, ++idx) // bug-testing, shorter and inaccurate loop
	{
		auto outer = tic();
		// const MatrixXd oldQ = Q;
		l1i.setZero();
		for (int k = 0; k < m; ++k) {
			// auto time = tic();
			l1ik.setZero();
			l1ik.head(X_nc) = X_uni.row(k).transpose() * yMinusLogi(i,k);
			l1ik.segment(X_nc, ZW_nc) = ZW.row(idx).transpose() * yMinusLogi(i,k);
			l1ik(ncov + k) = oneOverP(k);
			l1i += q(i,k) * l1ik;
			Q -= q(i,k) * l1ik * l1ik.transpose();
			// Rcout << "inner loop WALD "  << chrono::duration<double> (tic() - time).count() << endl;
		}
		Q += l1i * l1i.transpose();
		Rcout << "outer loop WALD "   << chrono::duration<double> (tic() - outer).count() << endl;
		// Rcout << Q.block(0,0, ncov, ncov) << endl << "Diff: " << (Q - oldQ).block(0,0, ncov, ncov).cwiseAbs().maxCoeff() << endl;

	}

	for (int k=0; k<m-1; ++k)
	{
		Q.block(0,ncov+k,ncov,1) -= Q.topRightCorner(ncov,1);
		for (int kk=k; kk<m-1; ++kk)
		{
			Q(ncov+k,ncov+kk) -= Q(ncov+k,dimQ-1) + Q(ncov+kk,dimQ-1) - Q(dimQ-1,dimQ-1);
		}
	}

	/**** augment the upper diagonal of Q **********************************************************************************************************/
	inv_profile_mat = Q.topLeftCorner(dimQ-1,dimQ-1).selfadjointView<Eigen::Upper>().ldlt().solve(MatrixXd::Identity(dimQ-1, dimQ-1));
	return inv_profile_mat.topLeftCorner(ncov,ncov);
} // WaldLogisticVarianceMLE0

// [[Rcpp::export]]
List TwoPhase_MLE0_logistic (const VectorXd& Y, const MatrixXd& X, const MatrixXd& ZW, const int& MAX_ITER, const double& TOL, const int& noSE)
{
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
	Create_Uni_Matrix_Row(X, X_index, X_uni, X_uni_ind, X_uni_n);
	/**** summarize observed distinct rows of X ****************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); // regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; // flag of none convergence in the estimation of regression coefficients
	bool flag_nonconvergence_cov; // flag of none convergence in the estimation of covariance matrix
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd LS_XtY(ncov);
	MatrixXd LS_XtX(ncov, ncov);
	VectorXd theta0(ncov);
	VectorXd p(m);
	VectorXd p0(m);
	MatrixXd q(n_minus_n2, m);
	VectorXd q_row_sum(n_minus_n2);
	RowVectorXd q_col_sum(m);
	MatrixXd P_theta(n_minus_n2, m);
	double tol;
	int iter, idx;
	VectorXd ZW_theta(n);
	VectorXd X_uni_theta(m);
	MatrixXd e_X_uni_theta(n_minus_n2, m);
	MatrixXd logi_X_uni_theta(n_minus_n2, m);
	MatrixXd logi2_X_uni_theta(n_minus_n2, m);
	VectorXd e_X_theta(n2);
	VectorXd logi_X_theta(n2);
	VectorXd logi2_X_theta(n2);
	string error_message;
	// /* RT test block */
	// time_t t1, t2;
	// /* RT test block */
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/

	/**** parameter initialization *****************************************************************************************************************/
	theta.setZero();
	theta0.setZero();
	for (int k=0; k<m; ++k)
	{
		p(k) = (double)X_uni_n(k) / (double)n2 ;
	}
	p0 = p;
	flag_nonconvergence = false;
	flag_nonconvergence_cov = false;
	/**** parameter initialization *****************************************************************************************************************/

	for (iter=0; iter<MAX_ITER; ++iter)
	{
		// /* RT test block */
		//	time(&t1);
		//	/* RT test block */

		/**** E-step *******************************************************************************************************************************/
		auto time = tic();
		/**** update P_theta ***********************************************************************************************************************/
		ZW_theta = ZW*theta.tail(ZW_nc);
		X_uni_theta = X_uni*theta.head(X_nc);
		for (int i=0; i<n_minus_n2; ++i)
		{
			for (int k=0; k<m; ++k)
			{
				e_X_uni_theta(i,k) = ZW_theta(i+n2) + X_uni_theta(k);
			}
		}
		e_X_uni_theta = e_X_uni_theta.array().exp();
        logi_X_uni_theta = e_X_uni_theta.array()/(e_X_uni_theta.array()+1.);
        logi2_X_uni_theta = logi_X_uni_theta.array()/(e_X_uni_theta.array()+1.);

		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			if (Y(idx) == 1.)
			{
				P_theta.row(i) = logi_X_uni_theta.row(i);
			}
			else if (Y(idx) == 0.)
			{
				P_theta.row(i) = 1.-logi_X_uni_theta.row(i).array();
			}
			else
			{
				std::ostringstream error_message;
				error_message << "Error: the " << i+n2 << "th value of Y is not 1 or 0!\n";
				stop(error_message.str());
			}
		}
		/**** update P_theta ***********************************************************************************************************************/
		Rcout << "update p theta "  << chrono::duration<double> (tic() - time).count() << endl;
		time = tic();
		/**** update q, q_row_sum ******************************************************************************************************************/
		for (int i=0; i<n_minus_n2; ++i)
		{
			for (int k=0; k<m; ++k)
			{
				q(i,k) = P_theta(i,k)*p(k);
			}
			// q.row(i) = P_theta.row(i).array() * p.array();
		}

		q_row_sum = q.rowwise().sum();
		for (int i=0; i<n_minus_n2; i++)
		{
			q.row(i) /= q_row_sum(i);
		}
		q_col_sum = q.colwise().sum();
		/**** update q, q_row_sum ******************************************************************************************************************/
		Rcout << "q, q_row_sum, "  << chrono::duration<double> (tic() - time).count() << endl;
		/**** E-step *******************************************************************************************************************************/


		/**** M-step *******************************************************************************************************************************/
		time = tic();
		/**** update theta *************************************************************************************************************************/
		for (int i=0; i<n2; i++)
		{
			e_X_theta(i) = ZW_theta(i)+X_uni_theta(X_uni_ind(i));
		}
		e_X_theta = e_X_theta.array().exp();
        logi_X_theta = e_X_theta.array()/(e_X_theta.array()+1.);
        logi2_X_theta = logi_X_theta.array()/(e_X_theta.array()+1.);

		LS_XtX.setZero();
		LS_XtY.setZero();
		const VectorXd yMinusLogi = Y - logi_X_theta;
		for (int i=0; i<n2; i++)
		{
			// LS_XtY.head(X_nc).noalias() += X.row(i).transpose()*(Y(i)-logi_X_theta(i));
			// LS_XtY.tail(ZW_nc).noalias() += ZW.row(i).transpose()*(Y(i)-logi_X_theta(i));
			LS_XtY.head(X_nc).noalias() += X.row(i).transpose() * yMinusLogi(i);
			LS_XtY.tail(ZW_nc).noalias() += ZW.row(i).transpose()* yMinusLogi(i);

			LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += (X.row(i).transpose())*X.row(i)*logi2_X_theta(i);
			LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += (X.row(i).transpose())*ZW.row(i)*logi2_X_theta(i);
			LS_XtX.bottomRightCorner(ZW_nc,ZW_nc).noalias() += (ZW.row(i).transpose())*ZW.row(i)*logi2_X_theta(i);
		}
		for (int i=0; i<n_minus_n2; ++i)
		{
			idx = i+n2;
			for (int k=0; k<m; ++k)
			{
				LS_XtY.head(X_nc).noalias() += q(i,k)*X_uni.row(k).transpose()*(Y(idx)-logi_X_uni_theta(i,k));
				LS_XtY.tail(ZW_nc).noalias() += q(i,k)*ZW.row(idx).transpose()*(Y(idx)-logi_X_uni_theta(i,k));

				LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += q(i,k)*(X_uni.row(k).transpose())*X_uni.row(k)*logi2_X_uni_theta(i,k);
				LS_XtX.topRightCorner(X_nc,ZW_nc).noalias() += q(i,k)*(X_uni.row(k).transpose())*ZW.row(idx)*logi2_X_uni_theta(i,k);
				LS_XtX.bottomRightCorner(ZW_nc,ZW_nc).noalias() += q(i,k)*(ZW.row(idx).transpose())*ZW.row(idx)*logi2_X_uni_theta(i,k);
			}
		}

		theta = LS_XtX.selfadjointView<Eigen::Upper>().ldlt().solve(LS_XtY);
		theta += theta0;
		/**** update theta *************************************************************************************************************************/
		Rcout << "update theta "  << chrono::duration<double> (tic() - time).count() << endl;
		time = tic();
		/**** update p *****************************************************************************************************************************/
		for (int k=0; k<m; ++k)
		{
			p(k) = X_uni_n(k)+0.;
		}
		p += q_col_sum.transpose();
		p /= n+0.;
		/**** update p *****************************************************************************************************************************/
		Rcout << "update p " << chrono::duration<double> (tic() - time).count() << endl;
		/**** M-step *******************************************************************************************************************************/


		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/
		tol = (theta-theta0).array().abs().sum();
		tol += (p-p0).array().abs().sum();
		/**** calculate the sum of absolute differences between estimates in the current and previous iterations ***********************************/

		/**** update parameters ********************************************************************************************************************/
		theta0 = theta;
		p0 = p;
		/**** update parameters ********************************************************************************************************************/

		/**** check convergence ********************************************************************************************************************/
		if (tol < TOL)
		{
			break;
		}
		/**** check convergence ********************************************************************************************************************/

		// /* RT test block */
		//	time(&t2);
		//	Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// /* RT test block */
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
		cov_theta.setConstant(-999.);
	}
	else if (noSE)
	{
		flag_nonconvergence_cov = true;
		cov_theta.setConstant(-999.);
	}
	else
	{
		cov_theta = WaldLogisticVarianceMLE0(LS_XtX, p, logi_X_uni_theta, Y, ZW, X_uni, q, n_minus_n2, X_nc, ZW_nc, m, ncov, n, n2);
	}
	/**** variance estimation **********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("theta") = theta,
						Named("cov_theta") = cov_theta,
						Named("flag_nonconvergence") = flag_nonconvergence,
						Named("flag_nonconvergence_cov") = flag_nonconvergence_cov);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // TwoPhase_MLE0_logistic
