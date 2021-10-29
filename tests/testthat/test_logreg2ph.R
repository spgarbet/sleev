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
	})