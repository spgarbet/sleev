library(testthat)
library(sleev)

test_that("MEXY", {
	skip_on_cran()
	skip_if(SKIP_CRAN_TESTS)

	args = commandArgs(TRUE)
	njob = 0
	rho = -.3
	p = 0.3
	hn_scale = 1
	nsieve = 20
	# wd = args[6]

	library(splines)

	NSIM = 20
	n = 1000
	n2 = 400
	alpha = 0.3
	beta = 0.4

	# setwd(wd)
	# fn_out = paste0(njob, ".Rdata")

	set.seed(12345)

	results = matrix(NA, nrow=NSIM, ncol=2)
	colnames(results) = c("Effect", "SE")

	nsim = 1
	while (nsim <= NSIM) {

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

	    # # histogram basis
	    # Bspline = matrix(NA, nrow=n, ncol=nsieve)
	    # cut_x_tilde = cut(simX_tilde, breaks=quantile(simX_tilde, probs=seq(0, 1, 1/nsieve)), include.lowest = TRUE)
	    # for (i in 1:nsieve) {
	    #     Bspline[,i] = as.numeric(cut_x_tilde == names(table(cut_x_tilde))[i])
	    # }
	    # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
	    # # histogram basis

	    # # linear basis
	    # Bspline = bs(simX_tilde, df=nsieve, degree=1, Boundary.knots=range(simX_tilde), intercept=TRUE)
	    # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
	    # # linear basis

	    # # quadratic basis
	    # Bspline = bs(simX_tilde, df=nsieve, degree=2, Boundary.knots=range(simX_tilde), intercept=TRUE)
	    # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
	    # # quadratic basis

	    # cubic basis
	    Bspline = bs(simX_tilde, df=nsieve, degree=3, Boundary.knots=range(simX_tilde), intercept=TRUE)
	    colnames(Bspline) = paste("bs", 1:nsieve, sep="")
	    # cubic basis

	    data = data.frame(Y_tilde=simY_tilde, X_tilde=simX_tilde, Y=simY, X=simX, Bspline)
	    ### generate data
	    expect_vector(data["Y"])
	    expect_vector(data["X"])
	    expect_vector(data["Y_tilde"])
	    expect_vector(data["X_tilde"])
	    expect_vector(data[colnames(Bspline)])

	    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1, verbose=TRUE)

	    if (sum(is.na(res$coefficients)) == 0) {
	        results[nsim,] = res$coefficients[2,1:2]
	        nsim = nsim+1
	        if (nsim%%10 == 0) {
	            print(paste(nsim, "replicates done."))
	        }
	    }
	}

	# Sep 22: 125 sec elapsed
	# $coefficients
	#            Estimate         SE Statistic      p-value
	# Intercept 0.2566482 0.03637235  7.056134 1.711964e-12
	# X         0.3811712 0.03559319 10.709107 0.000000e+00
	#
	# $sigma
	# [1] 0.9869545
	#
	# $covariance
	#               [,1]          [,2]
	# [1,]  1.322948e-03 -3.398867e-05
	# [2,] -3.398867e-05  1.266875e-03

	expect_equal(as.vector(res[["coefficients"]]), c(0.256648214727935,0.381171236418505,0.0363723532152385,0.0355931871009512,7.05613445490269,10.7091066427237,1.71196390397199e-12,0))
	expect_equal(res[["sigma"]], 0.986954535271583)
	expect_equal(as.vector( res[["covariance"]]), c(0.00132294807841407,-3.39886726125453e-05,-3.39886726125453e-05,0.00126687496800332))
	expect_true(res[["converge"]])
	expect_true(res[["converge_cov"]])
	})
