library(testthat)
library(sleev)


test_that("check use_default_values", {


	expect_equal(rho , -.3)
	expect_equal(p , 0.3)
	expect_equal(hn_scale , 1)
	expect_equal(nsieve , 20)
	expect_equal(NSIM , 20)
	expect_equal(n , 1000)
	expect_equal(n2 , 400)
	expect_equal(alpha , 0.3)
	expect_equal(beta , 0.4)
	})



test_that("MEXY loop", {
	skip_on_cran()
	skip_if(SKIP_CRAN_TESTS)


	# setwd(wd)
	# fn_out = paste0(njob, ".Rdata")

	set.seed(12345)

	results = matrix(NA, nrow=NSIM, ncol=2)
	colnames(results) = c("Effect", "SE")

	nsim = 1
	while (nsim <= NSIM) {

	    ### generate data
	    dat <- generate_data()
	    data <- dat$data
	    Bspline <- dat$Bspline

	    ### verify data
	    expect_vector(data["Y"])
	    expect_vector(data["X"])
	    expect_vector(data["Y_tilde"])
	    expect_vector(data["X_tilde"])
	    expect_vector(data[colnames(Bspline)])

	    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1, verbose=FALSE)

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

	## test all results
	expect_equal(as.vector(results), c(0.416922224791497,0.421737071865375,0.449369990588884,0.404940432955734,0.441158248372993,0.420228349283308,0.408144639683189,0.353816560047355,0.374311184655038,0.35189760767373,0.364188118169361,0.430859513855779,0.475390564039591,0.466674129567292,0.390630455867759,0.400706508262216,0.385259563151863,0.413198746826857,0.408625387823128,0.381171236418505,0.0374935035224319,0.0400514928365967,0.0375980846931493,0.0365689812070292,0.0344668558168398,0.0386203865132133,0.0393329246575616,0.0407566656906075,0.0368147916449737,0.0386479005173487,0.0379292995598047,0.03914357919284,0.0366241801208748,0.0405623820990424,0.0358255980316633,0.0403228646102295,0.0372884327822953,0.0373447475585087,0.0388191019382307,0.0355931871009512))

	## tests last res
	expect_equal(as.vector(res[["coefficients"]]), c(0.256648214727935,0.381171236418505,0.0363723532152385,0.0355931871009512,7.05613445490269,10.7091066427237,1.71196390397199e-12,0))
	expect_equal(res[["sigma"]], 0.986954535271583)
	expect_equal(as.vector( res[["covariance"]]), c(0.00132294807841407,-3.39886726125453e-05,-3.39886726125453e-05,0.00126687496800332))
	expect_true(res[["converge"]])
	expect_true(res[["converge_cov"]])
	})


test_that("single iteration cubic", {


	set.seed(12345)

    ### generate data
    dat <- generate_data()
    data <- dat$data
    Bspline <- dat$Bspline

    ### verify data
    expect_false(is.null(data))
    expect_false(is.null(Bspline))
    expect_vector(data["Y"])
    expect_vector(data["X"])
    expect_vector(data["Y_tilde"])
    expect_vector(data["X_tilde"])
    expect_vector(data[colnames(Bspline)])

    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1)

    expect_equal(as.vector(res[["coefficients"]]), c(0.288000325625334,0.416922224791497,0.036445693699495,0.0374935035213232,7.90217708572041,11.1198523913452,2.77555756156289e-15,0))
    expect_equal(res[["sigma"]], 1.01304717554374)
    expect_equal(as.vector( res[["covariance"]]), c(0.00132828858923741,-0.000111130454918251,-0.000111130454918251,0.00140576280630348))
    expect_true(res[["converge"]])
    expect_true(res[["converge_cov"]])
    })

test_that("single iteration histogram", {


	set.seed(12345)

    ### generate data
    dat <- generate_data("histogram")
    data <- dat$data
    Bspline <- dat$Bspline

    ### verify data
    expect_false(is.null(data))
    expect_false(is.null(Bspline))
    expect_vector(data["Y"])
    expect_vector(data["X"])
    expect_vector(data["Y_tilde"])
    expect_vector(data["X_tilde"])
    expect_vector(data[colnames(Bspline)])

    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1)

    expect_equal(as.vector(res[["coefficients"]]), c(0.282663936031731,0.414875913181448,0.0359350685141081,0.0374334492112151,7.86596346465173,11.0830265958273,3.66373598126302e-15,0))
    expect_equal(res[["sigma"]], 1.02778853831618)
    expect_equal(as.vector( res[["covariance"]]), c(0.00129132914911364,-8.59902916973147e-05,-8.59902916973147e-05,0.00140126311984862))
    expect_true(res[["converge"]])
    expect_true(res[["converge_cov"]])
    })

test_that("single iteration quadratic", {


	set.seed(12345)

    ### generate data
    dat <- generate_data("quadratic")
    data <- dat$data
    Bspline <- dat$Bspline

    ### verify data
    expect_false(is.null(data))
    expect_false(is.null(Bspline))
    expect_vector(data["Y"])
    expect_vector(data["X"])
    expect_vector(data["Y_tilde"])
    expect_vector(data["X_tilde"])
    expect_vector(data[colnames(Bspline)])

    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1)

    expect_equal(as.vector(res[["coefficients"]]), c(0.285986700184186,0.415752239593525,0.0363612977093676,0.0377064792460855,7.8651400857596,11.0260158971667,3.66373598126302e-15,0))
    expect_equal(res[["sigma"]], 1.01656868034353)
    expect_equal(as.vector( res[["covariance"]]), c(0.00132214397110926,-0.000108858715432434,-0.000108858715432434,0.00142177857713548))
    expect_true(res[["converge"]])
    expect_true(res[["converge_cov"]])
    })

test_that("single iteration linear", {


	set.seed(12345)

    ### generate data
    dat <- generate_data("linear")
    data <- dat$data
    Bspline <- dat$Bspline

    ### verify data
    expect_false(is.null(data))
    expect_false(is.null(Bspline))
    expect_vector(data["Y"])
    expect_vector(data["X"])
    expect_vector(data["Y_tilde"])
    expect_vector(data["X_tilde"])
    expect_vector(data[colnames(Bspline)])

    res = smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1)

    expect_equal(as.vector(res[["coefficients"]]), c(0.283549671920374,0.419193942728094,0.036188096643911,0.0377345448037876,7.83544033029666,11.1090234401348,4.66293670342566e-15,0))
    expect_equal(res[["sigma"]], 1.01921044143961)
    expect_equal(as.vector( res[["covariance"]]), c(0.00130957833870905,-0.000102885082209359,-0.000102885082209359,0.00142389587154905))
    expect_true(res[["converge"]])
    expect_true(res[["converge_cov"]])
    })

test_that("missing values", {

	set.seed(12345)

    ### generate data
    dat <- generate_data()
    data <- dat$data
    Bspline <- dat$Bspline


	# missing data
	expect_error(smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), hn_scale=0.1), "No dataset is provided!")

    # missing Y_unval
    expect_error(smle_MEXY(Y="Y", X="X", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1), "The error-prone response Y_unval is not specified!")

    # missing X_unval
    expect_error(smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1), "The error-prone covariates X_unval is not specified!")

    # missing Bspline
    expect_error(smle_MEXY(Y="Y", X="X", Y_unval="Y_tilde", X_unval="X_tilde", data=data, hn_scale=0.1), "The B-spline basis is not specified!")

    # missing Y
    expect_error(smle_MEXY(X="X", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1), "The accurately measured response Y is not specified!")

    # missing X
    expect_error(smle_MEXY(Y="Y", Y_unval="Y_tilde", X_unval="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1), "The validated covariates in the second-phase are not specified!")

    })



