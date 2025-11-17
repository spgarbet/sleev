# From main vignette
data("mock.vccc")

mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val / 10)
mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval / 10)
mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)
mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)

test_that("logistic2ph matches validated results", {
  sn <- 20
  expect_no_error(
    data.logistic <-
      spline2ph(x = "CD4_unval_sq10", size = 20, degree = 3,
                data = mock.vccc, group = "Prior_ART",
                split_group = TRUE)
  )

  expect_no_error({
    start.time <- Sys.time()
    res_logistic <- logistic2ph(y = "ADE_val", y_unval = "ADE_unval",
                                x = "CD4_val_sq10", x_unval = "CD4_unval_sq10",
                                z = "Prior_ART", data = data.logistic,
                                hn_scale = 1/2, se = TRUE, tol = 1e-04,
                                max_iter = 1000, verbose = FALSE)
    paste0("Run time: ", round(difftime(Sys.time(), start.time,
                                        units = "secs"), 3), " sec")
  })

  expect_no_error(res_logistic_coef <- summary(res_logistic)$coefficients)

  expect_close(
    res_logistic_coef,
    structure(c(-0.845963269605459, -0.541768855276144, -0.213722269892383,
                 0.304394570060743,  0.0628914007750632, 0.270579821947136,
                -2.77916675529608,  -8.61435504058544,  -0.789867730543996,
                 0.00544985400930287,0,                  0.429605018725176),
              dim = 3:4,
              dimnames = list(c("Intercept", "CD4_val_sq10", "Prior_ART"),
                              c("Estimate",  "SE", "Statistic", "p-value"))),
    1e-3)
})
