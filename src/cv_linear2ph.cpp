#include <RcppEigen.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "utility.h"

using namespace Rcpp;

//' @noRd
// [[Rcpp::export(.TwoPhase_MLE0_MEXY_CV_loglik)]]
List TwoPhase_MLE0_MEXY_CV_loglik (
	const Eigen::Map<Eigen::VectorXd>& Y_tilde, 
	const Eigen::Map<Eigen::MatrixXd>& X_tilde, 
	const Eigen::Map<Eigen::VectorXd>& Y, 
	const Eigen::Map<Eigen::MatrixXd>& X,
	const Eigen::Map<Eigen::MatrixXd>& Z, 
	const Eigen::Map<Eigen::MatrixXd>& Bspline, 
	const int& MAX_ITER, 
	const double& TOL, 
	const Eigen::Map<Eigen::VectorXd>& Train) 
{
	// /*#############################################################################################################################################*/
	// /**** pass arguments from R to cpp *************************************************************************************************************/
	// const MapVecd Y_tilde(as<MapVecd>(Y_tilde_R));
	// const MapMatd X_tilde(as<MapMatd>(X_tilde_R));
	// const MapVecd Y(as<MapVecd>(Y_R));
	// const MapMatd X(as<MapMatd>(X_R));
	// const MapMatd Z(as<MapMatd>(Z_R));
	// const MapMatd Bspline(as<MapMatd>(Bspline_R));
	// const int MAX_ITER = IntegerVector(MAX_ITER_R)[0];
	// const double TOL = NumericVector(TOL_R)[0];
	// const MapVecd Train(as<MapVecd>(Train_R));
	// /**** pass arguments from R to cpp *************************************************************************************************************/
	// /*#############################################################################################################################################*/
	
	
	
	/*#############################################################################################################################################*/
	/**** some useful constants ********************************************************************************************************************/
	const int n = Y_tilde.size();  // number of subjects in the first phase
	const int n2 = Y.size(); // number of subjects in the second phase
	const int n_minus_n2 = n-n2; // number of subjects not selected in the second phase
	const int Z_nc = Z.cols(); // number of inexpensive covariates
	const int X_nc = X.cols(); // number of expensive covariates X
	const int ncov = X_nc+Z_nc; // number of all covariates
	const int s = Bspline.cols(); // number of B-spline functions
	const int n_train = int(Train.sum());
	/**** some useful constants ********************************************************************************************************************/
	/*#############################################################################################################################################*/



	/*#############################################################################################################################################*/
	/**** error terms W+U **************************************************************************************************************************/
	const int WU_nc = 1+X_nc;
	MatrixXd WU(n2, WU_nc);
	WU.col(0) = Y_tilde.head(n2)-Y;
	WU.rightCols(X_nc) = X_tilde.topRows(n2)-X;
	/**** error terms W+U **************************************************************************************************************************/
	/*#############################################################################################################################################*/	
	

	
	/*#############################################################################################################################################*/	
	/**** summarize observed distinct rows of WU ***************************************************************************************************/
	// m: number of distinct rows of WU
	// WU_uni_ind(n2): row index of rows of WU in WU_uni
	// WU_uni(m, WU_nc): distinct rows of WU
	// WU_uni_n(m): count of appearances of each distinct row of WU	
	VectorXi WU_index = indexx_Matrix_Row(WU);
	int m = Num_Uni_Matrix_Row(WU, WU_index);	
	VectorXi WU_uni_ind(n2);		
	MatrixXd WU_uni(m, WU_nc); 
	VectorXi WU_uni_n(m);	
	Create_Uni_Matrix_Row(WU, WU_index, WU_uni, WU_uni_ind, WU_uni_n);
	/**** summarize observed distinct rows of WU ***************************************************************************************************/
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
	MatrixXd TrainDiag = Train.asDiagonal();
	
	// t(Y)*Y
	double LS_YtY_static = (TrainDiag.topLeftCorner(n2, n2)*Y).squaredNorm();
	LS_YtY_static += (TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*Y_tilde.tail(n_minus_n2)).squaredNorm(); 
	
	// t(X,Z)*Y
	VectorXd LS_XtY_static(ncov); 
	LS_XtY_static.head(X_nc) = X.transpose()*TrainDiag.topLeftCorner(n2, n2)*Y;
	LS_XtY_static.head(X_nc).noalias() += X_tilde.bottomRows(n_minus_n2).transpose()*TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*Y_tilde.tail(n_minus_n2);
	LS_XtY_static.tail(Z_nc) = Z.topRows(n2).transpose()*TrainDiag.topLeftCorner(n2, n2)*Y;
	LS_XtY_static.tail(Z_nc).noalias() += Z.bottomRows(n_minus_n2).transpose()*TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*Y_tilde.tail(n_minus_n2);
	
	// t(X,Z)*(X,Z)
	MatrixXd LS_XtX_static(ncov, ncov); 
	LS_XtX_static.topLeftCorner(X_nc,X_nc) = X.transpose()*TrainDiag.topLeftCorner(n2, n2)*X;
	LS_XtX_static.topLeftCorner(X_nc,X_nc).noalias() += X_tilde.bottomRows(n_minus_n2).transpose()*TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*X_tilde.bottomRows(n_minus_n2);	
	LS_XtX_static.topRightCorner(X_nc,Z_nc) = X.transpose()*TrainDiag.topLeftCorner(n2, n2)*Z.topRows(n2);
	LS_XtX_static.topRightCorner(X_nc,Z_nc).noalias() += X_tilde.bottomRows(n_minus_n2).transpose()*TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*Z.bottomRows(n_minus_n2);
	LS_XtX_static.bottomRightCorner(Z_nc,Z_nc) = Z.transpose()*TrainDiag*Z;	
	
	// p
	MatrixXd p_static(m, s); 
	p_static.setZero();
	for (int i=0; i<n2; i++) 
	{
		if (Train(i) == 1.) {
			p_static.row(WU_uni_ind(i)) += Bspline.row(i);
		}
	}	
	/**** some fixed quantities in the EM algorithm ************************************************************************************************/
	/*#############################################################################################################################################*/
	

	
	/*#############################################################################################################################################*/
	/**** output ***********************************************************************************************************************************/
	VectorXd theta(ncov); // regression coefficients
	MatrixXd cov_theta(ncov, ncov); // covariance matrix
	bool flag_nonconvergence; // flag of none convergence in the estimation of regression coefficients
	double sigma_sq; // residual variance
	double loglik = -999.;
	/**** output ***********************************************************************************************************************************/
	/*#############################################################################################################################################*/
	

		
	/*#############################################################################################################################################*/
	/**** temporary variables **********************************************************************************************************************/
	VectorXd LS_XtY(ncov);
	MatrixXd LS_XtX(ncov, ncov);
	VectorXd theta0(ncov);
	double LS_YtY, sigma_sq0;
	MatrixXd p(m, s); 
	MatrixXd p0(m, s);
	RowVectorXd p_col_sum(s);	
	MatrixXd q(n_minus_n2, m);
	VectorXd q_row_sum(n_minus_n2);
	RowVectorXd q_col_sum(m);
	MatrixXd pB(m_B, m);
	MatrixXd P_theta(n_minus_n2, m);
	double tol;
	int iter, idx;
	time_t t1, t2;
	/**** temporary variables **********************************************************************************************************************/
	/*#############################################################################################################################################*/
	

	
	/*#############################################################################################################################################*/
	/**** EM algorithm *****************************************************************************************************************************/
	
	/**** parameter initialization *****************************************************************************************************************/
	theta.setZero();
	theta0.setZero();
	sigma_sq = sigma_sq0 = Var(TrainDiag*Y_tilde)*n/n_train;

	p_col_sum = p_static.colwise().sum();
	for (int j=0; j<s; j++) 
	{
		if (p_col_sum(j) != 0.)
		{
			p.col(j) = p_static.col(j)/p_col_sum(j);
		} 
		else 
		{
			p.col(j).setZero();
		}
	}
	p0 = p;
	
	flag_nonconvergence = false;
	/**** parameter initialization *****************************************************************************************************************/
	
	for (iter=0; iter<MAX_ITER; iter++) 
	{
		time(&t1);
		Rcpp::checkUserInterrupt();
		
		/**** E-step *******************************************************************************************************************************/
		
		/**** update pB ****************************************************************************************************************************/
		pB = Bspline_uni*p.transpose();
		/**** update pB ****************************************************************************************************************************/
						
		/**** update P_theta ***********************************************************************************************************************/
		P_theta.col(0) = Y_tilde.tail(n_minus_n2);
		P_theta.col(0).noalias() -= X_tilde.bottomRows(n_minus_n2)*theta.head(X_nc);
		P_theta.col(0).noalias() -= Z.bottomRows(n_minus_n2)*theta.tail(Z_nc);
		for (int k=1; k<m; k++) 
		{
			P_theta.col(k) = P_theta.col(0);
		}		
		for (int k=0; k<m; k++) 
		{
			P_theta.col(k).noalias() += VectorXd::Constant(n_minus_n2, -WU_uni(k,0)+(WU_uni.block(k,1,1,X_nc)*theta.head(X_nc))(0,0));
		}		
		P_theta = P_theta.array().square();
		P_theta /= -2.*sigma_sq;
		P_theta = P_theta.array().exp();
		P_theta = (TrainDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*P_theta).eval();
		/**** update P_theta ***********************************************************************************************************************/
		
		/**** update q, q_row_sum ******************************************************************************************************************/
		q.setZero();
		for (int i=0; i<n_minus_n2; i++) 
		{
			if (Train(i+n2) == 1.)
			{
				for (int k=0; k<m; k++) 
				{
					q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
				}				
			}
		}
		q_row_sum = q.rowwise().sum();		
		for (int i=0; i<n_minus_n2; i++) 
		{
			if (Train(i+n2) == 1.)
			{
				q.row(i) /= q_row_sum(i);
			} 
			else
			{
				q.row(i).setZero();
				q_row_sum(i) = 0.;
			}
		}
		q_col_sum = q.colwise().sum();
		/**** update q, q_row_sum ******************************************************************************************************************/
		
		/**** E-step *******************************************************************************************************************************/
		
		
		/**** M-step *******************************************************************************************************************************/
		
		/**** update theta and sigma_sq ************************************************************************************************************/
		LS_YtY = LS_YtY_static;		
		LS_XtX = LS_XtX_static;
		LS_XtY = LS_XtY_static;
	
		for (int k=0; k<m; k++) 
		{
			LS_YtY += q_col_sum(k)*WU_uni(k,0)*WU_uni(k,0);
			LS_XtX.topLeftCorner(X_nc,X_nc).noalias() += q_col_sum(k)*WU_uni.block(k,1,1,X_nc).transpose()*WU_uni.block(k,1,1,X_nc);
			LS_XtY.head(X_nc).noalias() += q_col_sum(k)*WU_uni.block(k,1,1,X_nc).transpose()*WU_uni(k,0);		
		}
		
		for (int i=0; i<n_minus_n2; i++) 
		{
			idx = i+n2;
			if (Train(idx) == 1.) 
			{
				for (int k=0; k<m; k++) 
				{
					LS_YtY -= 2.*q(i,k)*Y_tilde(idx)*WU_uni(k,0);
					LS_XtX.topLeftCorner(X_nc,X_nc).noalias() -= q(i,k)*(WU_uni.block(k,1,1,X_nc).transpose()*X_tilde.row(idx)+X_tilde.row(idx).transpose()*WU_uni.block(k,1,1,X_nc));
					LS_XtX.topRightCorner(X_nc,Z_nc).noalias() -= q(i,k)*WU_uni.block(k,1,1,X_nc).transpose()*Z.row(idx);
					LS_XtY.head(X_nc).noalias() -= q(i,k)*(WU_uni.block(k,1,1,X_nc).transpose()*Y_tilde(idx)+X_tilde.row(idx).transpose()*WU_uni(k,0));
					LS_XtY.tail(Z_nc).noalias() -= q(i,k)*Z.row(idx).transpose()*WU_uni(k,0);
				}				
			}
		}
		
		theta = LS_XtX.selfadjointView<Eigen::Upper>().ldlt().solve(LS_XtY);		
		sigma_sq = LS_XtY.transpose()*theta;
		sigma_sq = (LS_YtY-sigma_sq)/n_train;
		/**** update theta and sigma_sq ************************************************************************************************************/
		
		/**** update p *****************************************************************************************************************************/
		p.setZero();		
		for (int i=0; i<n_minus_n2; i++) 
		{
			if (Train(i+n2) == 1.) 
			{
				for (int k=0; k<m; k++) 
				{
					for (int j=0; j<s; j++) 
					{
						p(k,j) += Bspline_uni(Bspline_uni_ind(i+n2),j)*P_theta(i,k)/q_row_sum(i);
					}
				}				
			}
		}		
		p = p.array()*p0.array();
		p += p_static;
		p_col_sum = p.colwise().sum();
		for (int j=0; j<s; j++) 
		{
			if (p_col_sum(j) != 0.)
			{
				p.col(j) /= p_col_sum(j);
			}
			else
			{
				p.col(j).setZero();
			}
		}
		/**** update p *****************************************************************************************************************************/
		
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
		
		time(&t2);
		// /* RT's test code */
		// Rcout << iter << '\t' << difftime(t2, t1) << '\t' << tol << endl;
		// /* RT's test code end */
	}
	/**** EM algorithm *****************************************************************************************************************************/
	/*#############################################################################################################################################*/
	
	
	
	/*#############################################################################################################################################*/
	/**** predicted log-likelihood *****************************************************************************************************************/
	if (iter == MAX_ITER) 
	{
		flag_nonconvergence = true;
	}
	else 
	{
		VectorXd resi_n(n2);
		MatrixXd logp(m, s);
		VectorXd Test = (1.-Train.array()).matrix();
		MatrixXd TestDiag = Test.asDiagonal();
		const int n_test = n-n_train;
		
		/**** calculate the likelihood *************************************************************************************************************/
		double tmp;
		
		logp.setZero();
		for (int k=0; k<m; k++)
		{
			for (int j=0; j<s; j++)
			{
				if (p(k,j) > 0.)
				{
					logp(k,j) = log(p(k,j));
				}
			}
		}
		pB = Bspline_uni*logp.transpose();
		
		loglik = 0.;
		for (int i=0; i<n2; i++) {
			if (Test(i) == 1.) {
				loglik += pB(Bspline_uni_ind(i),WU_uni_ind(i));
			}	
		}
		pB = Bspline_uni*p.transpose();
		
		P_theta.col(0) = Y_tilde.tail(n_minus_n2);
		P_theta.col(0).noalias() -= X_tilde.bottomRows(n_minus_n2)*theta.head(X_nc);
		P_theta.col(0).noalias() -= Z.bottomRows(n_minus_n2)*theta.tail(Z_nc);
		for (int k=1; k<m; k++) 
		{
			P_theta.col(k) = P_theta.col(0);
		}		
		for (int k=0; k<m; k++) 
		{
			P_theta.col(k).noalias() += VectorXd::Constant(n_minus_n2, -WU_uni(k,0)+(WU_uni.block(k,1,1,X_nc)*theta.head(X_nc))(0,0));
		}		
		P_theta = P_theta.array().square();
		P_theta /= -2.*sigma_sq;
		P_theta = P_theta.array().exp();
		P_theta = (TestDiag.bottomRightCorner(n_minus_n2, n_minus_n2)*P_theta).eval();
		
		q.setZero();
		for (int i=0; i<n_minus_n2; i++) 
		{
			if (Test(i+n2) == 1.)
			{
				for (int k=0; k<m; k++) 
				{
					q(i,k) = P_theta(i,k)*pB(Bspline_uni_ind(i+n2),k);
				}				
			}
		}
		q_row_sum = q.rowwise().sum();		
		for (int i=0; i<n_minus_n2; i++) 
		{
			if (Test(i+n2) == 1.)
			{
				loglik += log(q_row_sum(i));
			}
		}			
		loglik += -log(2.*M_PI*sigma_sq)*n_test/2.;
		resi_n = Y;
		resi_n.noalias() -= X*theta.head(X_nc);
		resi_n.noalias() -= Z.topRows(n2)*theta.tail(Z_nc);
		
		tmp = (TestDiag.topLeftCorner(n2,n2)*resi_n).squaredNorm();
		tmp /= 2.*sigma_sq;
		loglik -= tmp;
		/**** calculate the likelihood *************************************************************************************************************/
	}
	/**** predicted log-likelihood *****************************************************************************************************************/
	/*#############################################################################################################################################*/
	
	/*#############################################################################################################################################*/
	/**** return output to R ***********************************************************************************************************************/
	return List::create(Named("pred_loglike") = loglik,
						Named("flag_nonconvergence") = flag_nonconvergence);
	/**** return output to R ***********************************************************************************************************************/
	/*#############################################################################################################################################*/
} // TwoPhase_MLE0_MEXY_CV_loglik
