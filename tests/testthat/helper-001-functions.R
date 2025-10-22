expect_close <- function(x, y, tol=1e-8)
  expect_equal(as.numeric(x), as.numeric(y), tolerance=tol)
