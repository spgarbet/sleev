library(testthat)
library(CombinedReg)

test_that("0 vector", {
	w_t <- c(0,0,0)
	result <- lengthenWT(w_t, 3)

	expect_equal(as.vector(result), c(1,1,1,0,0,0))

	})
test_that("1 vector", {
  w_t <- c(1,1)
  result <- lengthenWT(w_t, 3)

  expect_equal(as.vector(result), c(1,1,1,1,1))

})
test_that("NA", {
	w_t <- NA
	result <- lengthenWT(w_t, 3)

		expect_equal(as.vector(result), c(1,1,1, NA))
	})

test_that("regular vector", {
	w_t <- seq(1, 20)
	result <- lengthenWT(w_t, 40)

	expect_equal(as.vector(result), c(rep(1, 40), seq(1,20)))
	})

# test_that("string", {
#   w_t <- "str"
#   expect_error(lengthenWT(w_t, 5))
#
#   expect_warning(lengthenWT(w_t, 5))
#
# })

test_that("don't modify", {
  w_t <- c(1,2,3)
  result <- lengthenWT(w_t, 4, FALSE)

  expect_equal(result, c(1,2,3))
})

test_that("n = 0",{
  w_t <- c(1,2,3)
  result <- lengthen(w_t, 0)

  expect_equal(result, c(1,2,3))
})
