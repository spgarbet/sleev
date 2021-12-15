library(testthat)
library(sleev)

test_that("documented example", {
	skip_if(!require(MASS))
	skip_if(!require(splines))
	
	rho = 0.3
	p = 0.3
	n = 100
	n2 = 40
	alpha = 0.3
	beta = 0.4
	 
	### generate data
	simX = rnorm(n)
	epsilon = rnorm(n)
	simY = alpha+beta*simX+epsilon

	error = MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1, rho, rho, 1), nrow=2))
	 
	simS = rbinom(n, 1, p)
	simU = simS*error[,2]
	simW = simS*error[,1]
	simY_tilde = simY+simW
	simX_tilde = simX+simU
	 
	id_phase2 = sample(n, n2)
	 
	simY[-id_phase2] = NA
	simX[-id_phase2] = NA
	 
	# cubic basis
	nsieves = c(5, 10, 15, 20, 25, 30, 40, 50, 80)
	pred_loglike = rep(NA, length(nsieves))
	for (i in 1:length(nsieves)) {
	    nsieve = nsieves[i]
	    Bspline = splines::bs(simX_tilde, df=nsieve, degree=3, Boundary.knots=range(simX_tilde), intercept=TRUE)
	    colnames(Bspline) = paste("bs", 1:nsieve, sep="")
	    # cubic basis
	   
	    data = data.frame(Y_tilde=simY_tilde, X_tilde=simX_tilde, Y=simY, X=simX, Bspline)
	    ### generate data
	   
	    res = cv_linear2ph(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, nfolds = 5)
	    pred_loglike[i] = res$avg_pred_loglik
	  }
	 
	data.frame(nsieves, pred_loglike)

	})
