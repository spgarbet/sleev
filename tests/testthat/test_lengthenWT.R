library(testthat)
library(CombinedReg)

test_that("0 vector", {
	w_t <- c(0,0,0)
	result <- lengthenWT(w_t, 3)

	expect_equal(result, c(1,1,1,0,0,0))
	})