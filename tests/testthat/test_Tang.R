library(testthat)

# test_that("fit_Tang", {
#     ## (5) MLE -------------------------------------------------
#     ### Script: two-phase log-likelihood specification adapted from Tang et al. (2015)
#     ### Get the script here https://github.com/sarahlotspeich/logreg2ph/blob/master/simulations/Tang_twophase_loglik_binaryX.R
#     source("simulations/Tang_twophase_loglik_binaryX.R")
#     fit_Tang <- nlm(f = Tang_twophase_loglik,
#                     p = rep(0, 12),
#                     hessian = TRUE,
#                     Val = "V",
#                     Y_unval = "Ystar",
#                     Y_val="Y",
#                     X_unval = "Xbstar",
#                     X_val = "Xb",
#                     C = "Xa",
#                     data = sdat)
#     beta_mle <- fit_Tang$estimate[10]
#     se_mle <- sqrt(diag(solve(fit_Tang$hessian)))[10]
    
# })
