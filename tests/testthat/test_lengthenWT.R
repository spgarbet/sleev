library(testthat)
library(CombinedReg)

test_that("0 vector", {
	w_t <- c(0,0,0)
	result <- as.vector(lengthenWT(w_t, 3))

	expect_equal(result, c(1,1,1,0,0,0))

	})

test_that("empty vector", {
	w_t <- c()
	result <- as.vector( lengthenWT(w_t, 3))

	expect_equal(result, c(1,1,1))
	})

test_that("regular vector", {
	w_t <- seq(1, 20)
	result <- as.vector(lengthenWT(w_t, 40))

	expect_equal(result, c(rep(1, 40), seq(1,20)))
	})
