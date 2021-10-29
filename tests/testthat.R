# initialize all
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(tictoc)
# Rcpp::sourceCpp('src/coxph.cpp')
# Rcpp::sourceCpp('src/linear.cpp')
Rcpp::sourceCpp('src/linearMEXY.cpp')
Rcpp::sourceCpp('src/lmm.cpp')
# Rcpp::sourceCpp('src/logistic.cpp')
# Rcpp::sourceCpp('src/utility.cpp')
source('R/smle.R')
source('R/smle_lmm.R')
source('R/smle_MEXY.R')

### SMLE TESTING ###

# library(TwoPhaseReg)
# n = 2000
# n2 = 600
# true_beta = 0.3
# true_gamma = 0.4
# true_eta = 0.5
# seed = 12345
# r = 0.3
# N_SIEVE = 8
#
# #### Sieve with histogram bases
# set.seed(12345)
# U2 = runif(n)
# simX = runif(n)
# simZ = r*simX+U2
# simY = true_beta*simX+true_gamma*simZ+rnorm(n)
# order.simY = order(simY)
# phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
# Bspline_Z = matrix(NA, nrow=n, ncol=N_SIEVE)
# cut_z = cut(simZ, breaks=quantile(simZ, probs=seq(0, 1, 1/N_SIEVE)), include.lowest = TRUE)
# for (i in 1:N_SIEVE) {
#     Bspline_Z[,i] = as.numeric(cut_z == names(table(cut_z))[i])
# }
# colnames(Bspline_Z) = paste("bs", 1:N_SIEVE, sep="")
# dat = data.frame(Y=simY, X=simX, Z=simZ, Bspline_Z)
# dat[-phase2.id,"X"] = NA
# #
# # res = smle(Y="Y", X="X", Z="Z", Bspline_Z=colnames(Bspline_Z), data=dat)
# # res
#
# #### Sieve with linear bases
# library(splines)
# set.seed(12345)
# U1 = runif(n)
# U2 = runif(n)
# simX = runif(n)
# simZ_1 = r*simX+U1
# simZ_2 = r*simX+U2
# simY = true_beta*simX+true_gamma*simZ_1+true_eta*simZ_2+rnorm(n)
# order.simY = order(simY)
# phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
# bs1 = bs(simZ_1, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_1), intercept=TRUE)
# bs2 = bs(simZ_2, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_2), intercept=TRUE)
# Bspline_Z = matrix(NA, ncol=N_SIEVE^2, nrow=n)
# for (i in 1:ncol(bs1)) {
#     for (j in 1:ncol(bs2)) {
#         idx = i*N_SIEVE+j-N_SIEVE
#         Bspline_Z[,idx] = bs1[,i]*bs2[,j]
#     }
# }
# colnames(Bspline_Z) = paste("bs", 1:ncol(Bspline_Z), sep="")
# dat = data.frame(Y=simY, X=simX, Z1=simZ_1, Z2=simZ_2, Bspline_Z)
# dat[-phase2.id,"X"] = NA
# sumnum <- 20
# dat["Delta"] = c(rep(1, sumnum), rep(0, 2000-sumnum))



### SMLE_MEXY TESTING ###
# args = commandArgs(TRUE)
# njob = 0
# rho = -.3
# p = 0.3
# hn_scale = 1
# nsieve = 20
# # wd = args[6]
#
# library(TwoPhaseReg)
# library(MASS)
# library(splines)
#
# NSIM = 20
# n = 1000
# n2 = 400
# alpha = 0.3
# beta = 0.4
#
# # setwd(wd)
# # fn_out = paste0(njob, ".Rdata")
#
# set.seed(12345)
#
# results = matrix(NA, nrow=NSIM, ncol=2)
# colnames(results) = c("Effect", "SE")
#
# nsim = 1
# tic(NSIM)
# while (nsim <= NSIM) {
#
#     ### generate data
#     simX = rnorm(n)
#     epsilon = rnorm(n)
#     simY = alpha+beta*simX+epsilon
#     error = mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1, rho, rho, 1), nrow=2))
#
#     simS = rbinom(n, 1, p)
#     simU = simS*error[,2]
#     simW = simS*error[,1]
#     simY_tilde = simY+simW
#     simX_tilde = simX+simU
#
#     id_phase2 = sample(n, n2)
#
#     simY[-id_phase2] = NA
#     simX[-id_phase2] = NA
#
#     # # histogram basis
#     # Bspline = matrix(NA, nrow=n, ncol=nsieve)
#     # cut_x_tilde = cut(simX_tilde, breaks=quantile(simX_tilde, probs=seq(0, 1, 1/nsieve)), include.lowest = TRUE)
#     # for (i in 1:nsieve) {
#     #     Bspline[,i] = as.numeric(cut_x_tilde == names(table(cut_x_tilde))[i])
#     # }
#     # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
#     # # histogram basis
#
#     # # linear basis
#     # Bspline = bs(simX_tilde, df=nsieve, degree=1, Boundary.knots=range(simX_tilde), intercept=TRUE)
#     # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
#     # # linear basis
#
#     # # quadratic basis
#     # Bspline = bs(simX_tilde, df=nsieve, degree=2, Boundary.knots=range(simX_tilde), intercept=TRUE)
#     # colnames(Bspline) = paste("bs", 1:nsieve, sep="")
#     # # quadratic basis
#
#     # cubic basis
#     Bspline = bs(simX_tilde, df=nsieve, degree=3, Boundary.knots=range(simX_tilde), intercept=TRUE)
#     colnames(Bspline) = paste("bs", 1:nsieve, sep="")
#     # cubic basis
#
#     data = data.frame(Y_tilde=simY_tilde, X_tilde=simX_tilde, Y=simY, X=simX, Bspline)
#     ### generate data
#
#     tic("smle_mexy")
#     res = smle_MEXY(Y="Y", X="X", Y_tilde="Y_tilde", X_tilde="X_tilde", Bspline=colnames(Bspline), data=data, hn_scale=0.1)
#     toc()
#     if (sum(is.na(res$coefficients)) == 0) {
#         results[nsim,] = res$coefficients[2,1:2]
#         nsim = nsim+1
#         if (nsim%%10 == 0) {
#             print(paste(nsim, "replicates done."))
#         }
#     }
# }
# res
# toc()
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


### SMLE_LMM TESTING ###

r = 1
design = "ods"
goal = "both"
N_SIEVE = 20

# library(devtools)
# install_github("dragontaoran/TwoPhaseReg")
library(lme4)
library(splines)

NSIM = 2
nid = 1000
nid2 = 400
n_per_id = 5
true_beta1 = 0.3
true_beta2 = 0.3
true_alpha2 = 0.4
true_gamma1 = 0.5
true_gamma2 = 0.5
psi11 = 1
psi22 = 1

set.seed(12345)

results_est = matrix(NA, nrow=NSIM, ncol=6)
colnames(results_est) = c("intercept", "X", "Z", "Time", "Time_X", "Time_Z")
results_se = matrix(NA, nrow=NSIM, ncol=6)
colnames(results_se) = c("intercept", "X", "Z", "Time", "Time_X", "Time_Z")
results_vc = matrix(NA, nrow=NSIM, ncol=4)
tic(NSIM)
for (nsim in 1:NSIM) {

    ### generate data
    U2 = runif(nid)
    simX = runif(nid)
    simZ = r*simX+U2
    simT = rep(seq(-1, 1, length.out=n_per_id), nid)
    simb0 = rnorm(nid, sd=sqrt(psi11))
    simb1 = rnorm(nid, sd=sqrt(psi22))
    simY = rep(true_beta1*simX+true_gamma1*simZ, each=n_per_id)
    simY = simY+true_alpha2*simT
    simY = simY+true_beta2*rep(simX, each=n_per_id)*simT
    simY = simY+true_gamma2*rep(simZ, each=n_per_id)*simT
    simY = simY+rep(simb0, each=n_per_id)
    simY = simY+rep(simb1, each=n_per_id)*simT
    simY = simY+rnorm(nid*n_per_id)
    simID = rep(1:nid, each=n_per_id)
    # Bspline_Z = matrix(NA, nrow=nid, ncol=N_SIEVE)
    # cut_z = cut(simZ, breaks=quantile(simZ, probs=seq(0, 1, 1/N_SIEVE)), include.lowest = TRUE)
    # for (i in 1:N_SIEVE) {
    #     Bspline_Z[,i] = as.numeric(cut_z == names(table(cut_z))[i])
    # }
    # colnames(Bspline_Z) = paste("bs", 1:N_SIEVE, sep="")

    # linear basis
    Bspline_Z = bs(simZ, df=N_SIEVE, degree=1, Boundary.knots=range(simZ), intercept=TRUE)
    colnames(Bspline_Z) = paste("bs", 1:N_SIEVE, sep="")
    # linear basis

    data = data.frame(Y=simY, Time=simT, ID=simID, Z=rep(simZ, each=n_per_id), Bspline_Z[simID,])

    if (design == "ods") {

        Q = matrix(NA, nid, 2)
        for (i in 1:nid) {
            idxx = which(simID == i)
            tmp = lm(simY[idxx]~simT[idxx])
            Q[i,] = coef(tmp)
        }
        order_Q = matrix(NA, nid, 2)
        order_Q[,1] = order(Q[,1])
        order_Q[,2] = order(Q[,2])

    } else if (design == "blup") {

        res0 = lmer(Y~Z+Time+Z*Time+(Time|ID), data=data, REML=FALSE)
        Q = ranef(res0)$ID
        order_Q = matrix(NA, nid, 2)
        order_Q[,1] = order(Q[,1])
        order_Q[,2] = order(Q[,2])
    }

    if (goal == "intercept") {
        id_phase2 = c(order_Q[1:(nid2/2),1], order_Q[(nid-(nid2/2)+1):nid,1])
    } else if (goal == "slope") {
        id_phase2 = c(order_Q[1:(nid2/2),2], order_Q[(nid-(nid2/2)+1):nid,2])
    } else if (goal == "both") {
        id_phase2_1 = c(order_Q[1:(nid2/4),1], order_Q[(nid-(nid2/4)+1):nid,1])
        ids_remain = (1:nid)[-id_phase2_1]
        nid_remain = length(ids_remain)
        order_Q2 = order(Q[ids_remain,2])
        id_phase2_2 = ids_remain[c(order_Q2[1:(nid2/4)], order_Q2[(nid_remain-(nid2/4)+1):nid_remain])]
        id_phase2 = c(id_phase2_1, id_phase2_2)
    }

    simX[-id_phase2] = NA
    data$X = rep(simX, each=n_per_id)

    tic("smle_lmm")
    res = smle_lmm(Y="Y", Time="Time", ID="ID", X="X", Z="Z", Bspline_Z=colnames(Bspline_Z), data=data) # , ZT=TRUE
    toc()
    # results_est[nsim,] = t(res$coefficients[,1])
    # results_se[nsim,] = t(res$coefficients[,2])
    # results_vc[nsim,] = t(res$vc)

    if (nsim%%50 == 0) {
        print(paste(nsim, "replicates done."))
    }
}
res
toc()
# 2: 286.84 sec elapsed
# Sep 24: 219 sec
# $coefficients
#              Estimate         SE Statistic      p-value
# Intercept -0.03852952 0.08944463 -0.430764 6.666400e-01
# X          0.22212501 0.19811734  1.121179 2.622117e-01
# Z          0.61596572 0.13162632  4.679655 2.873584e-06
# Time       0.56303492 0.07585907  7.422117 1.152411e-13
# Time:X     1.04636764 0.13535209  7.730709 1.065814e-14
#
# $vc
# [1]  0.95772860 -0.00672124  0.96407605  0.97952556


# tic("smle")


# TwoPhase_GeneralSpline
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), data=dat,model="linear")
# 177.35 sec elapsed
# $coefficients
#            Estimate         SE Statistic      p-value
# Intercept 0.0538190 0.07287324 0.7385291 4.601930e-01
# X         0.3367835 0.09518781 3.5380951 4.030249e-04
# Z1        0.3808042 0.07657005 4.9732790 6.582983e-07
# Z2        0.4135827 0.07696735 5.3734818 7.723062e-08
# Aug 24 : 133.49 sec elapsed
# Sep 13 : 125.75 sec elapsed



# TwoPhase_MLE0
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=NULL, data=dat)
# $coefficients
# Estimate         SE  Statistic      p-value
# Intercept -0.06193395 0.08597975 -0.7203319 4.713207e-01
# X          0.46180716 0.13255814  3.4838084 4.943332e-04
# Z1         0.42190671 0.07844466  5.3783993 7.515098e-08
# Z2         0.44886314 0.07999735  5.6109751 2.011897e-08
#  1724.11 sec elapsed



# TwoPhase_GeneralSpline_logistic
# dat["Y"] = sample(c(0,1), nrow(dat), TRUE)
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), data=dat, model="logistic")
# $coefficients
#              Estimate        SE  Statistic   p-value
# Intercept -0.03133511 0.1527562 -0.2051314 0.8374694
# X          0.35195573 0.3071224  1.1459786 0.2518040
# Z1        -0.09435242 0.1658927 -0.5687557 0.5695219
# Z2        -0.10004598 0.1647833 -0.6071365 0.5437603
# Sep 2: 148.69 sec elapsed
# Sep 8: 116.66 sec elapsed


# TwoPhase_MLE0_logistic
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=NULL, data=dat, model="logistic")
# $coefficients
#              Estimate        SE   Statistic   p-value
# Intercept -0.01363803 0.1830567 -0.07450169 0.9406112
# X          0.13783433 0.2940132  0.46880314 0.6392104
# Z1        -0.02346967 0.1507358 -0.15570075 0.8762689jk
# Z2        -0.03550027 0.1522553 -0.23316271 0.8156351
# 2539.78 sec elapsed

# 100 iterations - abbr
#              Estimate        SE   Statistic   p-value
# Intercept -0.01363803 0.1535751 -0.08880363 0.9292380
# X          0.13783433 0.1623791  0.84884282 0.3959688
# Z1        -0.02346967 0.1495529 -0.15693225 0.8752982
# Z2        -0.03550027 0.1510185 -0.23507233 0.8141526
# 139.8 sec


# TwoPhase_GeneralSpline_coxph
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), Delta="Delta", data=dat, model="coxph")
# 20 Delta
# Estimate        SE  Statistic   p-value
# X   1.1585078 1.8492915  0.6264603 0.5310131
# Z1 -0.5853910 0.8287508 -0.7063534 0.4799684
# Z2 -0.2120672 0.8410371 -0.2521496 0.8009254
# 534.81 sec elapsed

# TwoPhase_MLE0_noZW_coxph
# res = smle(Y="Y", X="X",  Delta="Delta", data=dat, model="coxph")
# too long to measure


# TwoPhase_MLE0_coxph INCONSISTENT OUTPUT
# res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Delta="Delta", data=dat, model="coxph")
#      Estimate            SE      Statistic p-value
# X  -1.0218959 2.237591e-154 -4.566946e+153       0
# Z1 -0.2691035 2.237542e-154 -1.202675e+153       0
# Z2  0.1267448 3.731355e-159  3.396750e+157       0
# smle: 275.84 sec elapsed

# Estimate            SE      Statistic   p-value
# X  -1.0218959  7.026472e+58  -1.454351e-59 1.0000000
# Z1 -0.2691035 7.103049e+111 -3.788563e-113 1.0000000
# Z2  0.1267448  3.669300e+02   3.454195e-04 0.9997244
# 511 sec

# Estimate            SE      Statistic   p-value
# X  -1.0218959 2.237591e-154 -4.566946e+153 0.0000000
# Z1 -0.2691035 2.237542e-154 -1.202675e+153 0.0000000
# Z2  0.1267448  1.414214e+00   8.962211e-02 0.9285875

# Estimate            SE      Statistic  p-value
# X  -1.0218959 1.134131e+124 -9.010384e-125 1.000000
# Z1 -0.2691035 3.383034e+121 -7.954503e-123 1.000000
# Z2  0.1267448  1.433320e+00   8.842741e-02 0.929537

# Estimate            SE      Statistic   p-value
# X  -1.0218959 1.134131e+124 -9.010384e-125 1.0000000
# Z1 -0.2691035 3.383034e+121 -7.954503e-123 1.0000000
# Z2  0.1267448  1.414253e+00   8.961958e-02 0.9285895

# res
# toc()
# add log likelihood
# add aic calculation

