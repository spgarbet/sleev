library(testthat)
# library(CombinedReg)
library(logreg2ph)

test_that("logreg2ph Simulation 1", {

	## (6) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### Since Xb* and Xa are both binary, reduces to indicators --
	nsieve <- 4
	B <- matrix(0, nrow = N, ncol = nsieve)
	B[which(Xa == 0 & Xbstar == 0), 1] <- 1
	B[which(Xa == 0 & Xbstar == 1), 2] <- 1
	B[which(Xa == 1 & Xbstar == 0), 3] <- 1
	B[which(Xa == 1 & Xbstar == 1), 4] <- 1
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)
	smle <- logreg2ph(Y_unval = "Ystar",
		Y_val = "Y",
		X_unval = "Xbstar",
		X_val = "Xb",
		C = "Xa",
		Validated = "V",
		Bspline = colnames(B),
		data = sdat,
		noSE = FALSE,
		MAX_ITER = 1000,
		TOL = 1E-4)

	expect_equal(smle[["coeff"]]["coeff"], c(-0.587801556603343, -0.231437717735943, 0.139139424615871))
	expect_equal(smle[["coeff"]]["se"], c(0.159867158182108, 0.240522213900834, 0.203146591074255))

	expect_equal(smle[["outcome_err_coeff"]]["coeff"], c(-2.50954051127845, -0.763949870590926, 5.6328754487314, 0.50008005682218, -0.345208034572803))
	expect_equal(smle[["outcome_err_coeff"]]["se"], c(NA, NA, NA, NA, NA))

	expect_equal(smle[["Bspline_coeff"]], c(1,1,0.744983827508479,0.255016172491521,0.24073358078479,0.75926641921521,0.65244609543555,0.34755390456445,0.352216846922798,0.647783153077202))
	expect_equal(smle[["vcov"]], c(0.0255575082652231,-0.0291909056522233,-0.0106240625411377,-0.0291909056522233,0.0578509353797586,-0.000415041049293844,-0.0106240625411377,-0.000415041049293844,0.0412685374650904))
	expect_equal(smle[["od_loglik_at_conv"]], -847.635948985586)

	expect_true(smle[["converged"]])
	expect_true(smle[["se_converged"]])
})