---
title: 'sleev: An R Package for Semiparametric Likelihood Estimation with Errors in
  Variables'
tags:
  - R
  - measurement error
  - "two-phase sampling"
  - nonparametric statistics
  - "B-splines"
date: "16 September 2024"
output: pdf_document
bibliography: paper.bib
affiliations:
  - name: Department of Biostatistics, Vanderbilt University Medical Center, USA
    index: 1
  - name: Department of Statistical Sciences, Wake Forest University, USA
    index: 2
  - name: Brigham Young University, USA
    index: 3
  - name: Vanderbilt Genetics Institute, Vanderbilt University Medical Center, USA
    index: 4
authors:
  - name: Jiangmei Xiong
  - corresponding: yes
    affiliation: 1
  - name: Sarah C. Lotspeich
    affiliation: 2
  - name: Joey B. Sherrill
    affiliation: 3
  - name: Gustavo Amorim
    affiliation: 1
  - name: Bryan E. Shepherd
    affiliation: 1
  - name: Ran Tao
    affiliation: '1,4'
---

# Summary

Data with measurement error in the outcome, covariates, or both are not
uncommon, particularly with the increased use of routinely collected
data for biomedical research. In settings with error-prone data, often
only a subsample of study data is validated, also known as two-phase
studies. The sieve maximum likelihood estimator (SMLE), which combines
the error-prone data on all records with the validated data on a
subsample, is a highly efficient and robust estimator to analyze such
data. However, given their complexity, a computationally efficient and
user-friendly tool is needed to obtain SMLEs. R package sleev fills this
gap by making semiparametric likelihood-based inference using SMLEs for
error-prone two-phase data in settings with binary and continuous
outcomes. Functions from this package can be used to analyze data with
error-prone responses and covariates.

# Statement of Need

Routinely collected data are being used frequently in biomedical
research, such as electronic health records. However, these data tend to
be error-prone, and careless use of these data could lead to biased
estimates and misleading research findings (@duan2016empirical). To avoid invalid
study results, trained experts carefully verify and extract data
elements. However, it is usually only feasible to validate data for a
subset of records or variables. After validation, researchers have
error-prone pre-validation data for all records (phase one data)
error-free validated data on a subset of records (phase two data).
Analyses aims to combine the two types of data to obtain estimates that
have low bias and are as efficient as possible.

SMLE, which combines the error-prone data on all records with the
validated data on a subsample, is a highly efficient and robust
estimator to analyze such two-phase data (@tao2021efficient; @lotspeich2022efficient). Still, in
practice these estimators can be difficult to implement, as they involve
approximating nuisance conditional densities using B-splines (@schumaker2007spline)
and then maximizing the semiparametric likelihood via a sophisticated EM
algorithm (@tao2017efficient). R package sleev makes this method readily
applicable for practitioners in a user-friendly way. sleev integrates
and extends R packages logreg2ph and TwoPhaseReg, two primitive R
packages developed with the original methods papers (@tao2021efficient; @lotspeich2022efficient). These
two packages lacked proper documentation, were computationally slow, and
were difficult to use. To promote the use of SMLE in two-phase data,
extensive work has been done to create sleev, a computationally
efficient and user-friendly R package to analyze two-phase, error-prone
data. Specifically, in sleev we rewrote the core algorithms in C++ to
speed up the computation, and we unified the syntax across functions.

# SMLE for linear regression

In this section, we briefly introduce the SMLEs for linear regression.
The SMLEs for logistic regression are similar to linear regression and
described in **reference to github full paper**.Suppose that we want to fit
a standard linear regression model for a continuous outcome $Y$ and
covariates $\mathbf{X}$:
$Y = \alpha + \mathbf{\beta}^{T}\ \mathbf{X} + \epsilon$, where
$\epsilon\sim N({0,\sigma}^{2})$. Our goal is to obtain estimates of
${\mathbf{\theta} = (\alpha,\ \mathbf{\beta}^{T},\sigma^{2})}^{T}$. When
we have error-prone data, $Y$ and $\mathbf{X}$ are unobserved except for
a subset of validated records. For the majority of unvalidated records,
only the error-prone outcome $Y^{*} = Y + W$ and covariates
$\mathbf{X}^{\mathbf{*}} = \mathbf{X} + \mathbf{U}$ are observed in
place of $Y$ and $\mathbf{X}$, where $W$ and $\mathbf{U}$ are the errors
for the outcome and covariates, respectively. We assume that
$W,\ \mathbf{U}\mathbf{\bot}\epsilon$. With potential errors in our
data, a naive regression analysis using error-prone variables $Y^{*}$
and $\mathbf{X}^{\mathbf{*}}$ could render misleading results.

We assume that the joint density of the complete data
$\left( Y^{*},\mathbf{X}^{*},W,\mathbf{U} \right)$ takes the form

$$P\left( Y^{*},\mathbf{X}^{*},\ W,\ \mathbf{U} \right) = P\left( Y^{*}|\mathbf{X}^{*},\ W,\ \mathbf{U} \right)P\left( W,\ \mathbf{U|}\mathbf{X}^{*}\  \right)P\left( \mathbf{X}^{*} \right)$$

$$= P_{\mathbf{\theta}}\left( Y|X \right)P\left( W,\ \mathbf{U|}\mathbf{X}^{*} \right)P\left( \mathbf{X}^{*} \right)$$

where $P( \cdot )$ and $P\left( \cdot | \cdot \right)$ denote density
and conditional density functions, respectively.
$P_{\mathbf{\theta}}\left( Y|\mathbf{X} \right)$ then refers to the
conditional density function of
$Y = \alpha + \mathbf{X\beta} + \epsilon$. Denote the validation
indicator variable by $V$, with $V = 1$ indicating that a record was
validated and $V = 0$ otherwise. For records with $V = 0$, their
measurement errors $\left( W,\mathbf{U} \right)$ are missing, and
therefore their contributions to the log-likelihood can be obtained by
integrating out $W$ and $\mathbf{U}$. Let
$\left( Y_{i}^{*},\mathbf{X}_{i}^{\mathbf{*}},W_{i},\mathbf{U}_{i},V_{i},Y_{i},\mathbf{X}_{i} \right)$
for $i = 1,\ldots,n$ denote independent and identically distributed
realizations of
$\left( Y^{*},\mathbf{X}^{\mathbf{*}},W,\mathbf{U},V,Y,\mathbf{X} \right)$
in a sample of $n$ subjects. Then, the observed-data log-likelihood
takes the form

$$\sum_{i = 1}^{n}{V_{i}\{\log P_{\mathbf{\theta}}\left( Y_{i} \middle| \mathbf{X}_{i} \right) + \log{P(W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{*})\}}}$$

$$+ \sum_{i = 1}^{n}{\left( 1 - V_{i} \right)\log\left\{ \int\int P_{\mathbf{\theta}}\left( Y_{i}^{*} - w|\mathbf{X}_{i}^{*} - \mathbf{u} \right)P\left( w,\mathbf{u} \middle| \mathbf{X}_{i}^{\mathbf{*}} \right)dwd\mathbf{u} \right\}}.\ \ \ \ \ \ \ \ (1)$$

We estimate the unknown measurement error model,
$P\left( W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{\mathbf{*}} \right)$ using
B-spline sieves. Specifically, we approximate
$P\left( w_{i},\mathbf{u}_{i}|\mathbf{X}_{i}^{\mathbf{*}} \right)$ and
$\log P\left( W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{\mathbf{*}} \right)$
by
$\sum_{k = 1}^{m}I\left( w_{i} = w_{k},\mathbf{u}_{i} = \mathbf{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)p_{kj}$
and
$\sum_{k = 1}^{m}I\left( W_{i} = w_{k},\mathbf{U}_{i} = \mathbf{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)\log p_{kj}$,
respectively.
$\left\{ \left( w_{1},\mathbf{u}_{1} \right),...,\left( w_{m},\mathbf{u}_{m} \right) \right\}$
are the $m$ distinct observed $\left( W,\mathbf{U} \right)$ values,
$B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)$ is the $j$th
B-spline basis function of order $q$ evaluated at
$\mathbf{X}_{i}^{\mathbf{*}}$, $s_{n}$ is the dimension of the B-spline
basis, and $p_{kj}$ is the coefficient associated with
$B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)$ and
$\left( w_{k},\mathbf{u}_{k} \right)$. The log-likelihood in expression
(1) is now approximated by

$$\sum_{i = 1}^{n}{V_{i}\{\log P_{\mathbf{\theta}}\left( Y_{i} \middle| \mathbf{X}_{i} \right) + \sum_{k = 1}^{m}{I\left( W_{i} = w_{i},\ \mathbf{U}_{i} = \mathbf{u}_{k} \right)\sum_{j = 1}^{S_{n}}{B_{j}^{q}(\mathbf{X}_{i}^{*})}}\log{p_{kj}\}}}$$

$$+ \sum_{i = 1}^{n}{\left( 1 - V_{i} \right)\log\left\{ \sum_{k = 1}^{m}{P_{\mathbf{\theta}}\left( Y_{i}^{*} - w|\mathbf{X}_{i}^{*} - \mathbf{u} \right)}\sum_{j = 1}^{S_{n}}{B_{j}^{q}(\mathbf{X}_{i}^{*})}\log{p_{kj}\}} \right\}}.\ \ \ \ \ \ \ (2)$$

The maximization of expression (2) is carried out through an EM
algorithm to find the SMLEs $\widehat{\mathbf{\theta}}$ and
${\widehat{p}}_{kj}$. The covariance matrix of the SMLE
$\widehat{\mathbf{\theta}}$ is obtained through the method of profile
likelihood (@murphy2000profile). Full details on the SMLE method for logistic regression
with error-prone data, including its theoretical properties, can be
found in (@lotspeich2022efficient).

# Functionalities of the `sleev` R package

The sleev package provides a user-friendly way to obtain the SMLEs and
their standard errors. The sleev package includes two main functions:
linear2ph() and logistic2ph(), to fit linear and logistic regressions,
respectively, under two-phase sampling with an error-prone outcome and
covariates. The input arguments are similar for the two functions and
listed in Table 1. From Table 1, we see that in addition to the
arguments for error-prone and error-free outcome and covariates, the
user needs to specify the B-spline matrix
$B_{j}^{q}\left( \mathbf{X}_{i}^{*} \right)$ to be used in the
estimation of the error densities.

The sleev package included a dataset constructed to mimic data from the
Vanderbilt Comprehensive Care Clinic (VCCC) patient records from (@giganti2020accounting).
Table 1 lists the variables in this dataset to be used in subsequent
analyses.

# Example: Case study with mock data

We now illustrate how to obtain SMLEs using the sleev package with
dataset mock.vccc. Specifically, we show how to fit a linear regression
model in the presence of errors in both the outcome and covariates using
the linear2ph() function. Situations with more covariates and examples
with logistic regression are included in **reference to github full paper**.

This example fits a linear regression model with CD4 count at antiretroviral therapy (ART) 
initiation regressed on viral load (VL) at ART initiation, adjusting for sex at birth. 
Both CD4 and VL are error-prone, partially validated variables, whereas sex is error-free.

Because of skewness, we often transform both CD4 and VL. In our
analysis, CD4 was divided by 10 and square root transformed and VL
was $\log_{10}$ transformed:

```
mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val/10)

mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval/10)

mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)

mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)
```

To obtain the SMLEs, we first need to set up the B-spline basis for the covariates `VL_unval_l10` (the transformed error-prone VL
variable from phase one) and `Sex`. Here, we use a cubic B-spline basis with the `degree=3` argument in our
call to the `bs()` function from the `splines` package. The size
of the basis $s_{n}$ is set to be 20; this was specified through the
`df = 20` argument. The B-spline basis is set up separately for the two
`Sex` groups, and the size of the B-spline basis is assigned in proportion to the relative size of the two `Sex` groups. This allows the errors in `VL_unval_l10` to be
heterogeneous between males and females. The described B-spline basis is constructed as follows.  

```
n <- nrow(mock.vccc) # get the number of rows of the dataset
sn <- 20 # set the size of the B-spline basis
# portion the size of the B-spline basis according to the size of two sex groups
sex_ratio <- sum(mock.vccc$Sex==0)/n
sn_0 <- round(sn*sex_ratio, digits = 0)
sn_1 <- sn-sn_0
# create B-spline basis for each sex group separately through bs()
Bspline_0 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$Sex==0], df=sn_0, degree=3,
                          intercept=TRUE)
Bspline_1 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$Sex==1], df=sn_1, degree=3,
                          intercept=TRUE)
# create B-spline basis for this analysis by combining the two bases created above
Bspline <- matrix(data=0, nrow=n, ncol=sn)
Bspline[mock.vccc$Sex==0,1:(sn_0)] <- Bspline_0
Bspline[mock.vccc$Sex==1,(sn_0+1):sn] <- Bspline_1
# name B-spline basis matrix columns properly
colnames(Bspline) <- paste0("bs", 1:sn) 
```
Alternatively, if the
investigator has prior knowledge that the errors in `VL_unval_l10` are
likely to be homogeneous, one may fit a simpler model by not stratifying
the B-spline basis by `Sex`. As a final step before running the analysis, we add the B-spline basis
to the analysis dataset to create a single data frame:

```
data <- data.frame(cbind(mock.vccc, Bspline))
```
The SMLEs can be obtained by running function `linear2ph()`, as shown
in the code below. The fitted SMLEs are stored in a list object. Here,
we assign the fitted SMLEs to the variable name res_linear.This list
contains five slots: coefficients, covariance, sigma, converge, and
converge_cov. We should first check if the EM algorithms for estimating
the regression coefficients and their covariance matrix converged by
checking if `res_linear$converge` and `res_linear$converge_conv` are TRUE.

The `res_linear$coefficient` slot gives the regression coefficient
estimates and their corresponding standard error estimates,
$z$-statistics, and $p$-values.

```
res_linear <- linear2ph(Y_unval="CD4_unval_sq10", Y="CD4_val_sq10", X_unval="VL_unval_l10",
                        X="VL_val_l10", Z="Sex", Bspline=colnames(Bspline), data=data,
                        hn_scale = 1, noSE = FALSE, TOL = 1e-04, MAX_ITER = 1000,
                        verbose = FALSE)
```

Check convergence:
```
c(res_linear$converge, res_linear$converge_cov)
```
Result coefficients:
```
res_linear$coefficients

# Estimate SE Statistic p-value
# Intercept 4.821 0.1587 30.39 0.000000
# VL_val_l10 -0.141 0.0398 -3.55 0.000389
# Sex 0.273 0.1089 2.51 0.012229
```

# Acknowledgement

This research was supported by the National Institute of Health grants
R01AI131771, R01HL094786, and P30AI110527 and the 2022 Biostatistics
Faculty Development Award from the Department of Biostatistics at
Vanderbilt University Medical Center.

# Reference


Table 1: Main arguments for linear2ph() and logistic2ph()

-----------------------------------------------------------------------
Argument            Description
------------------- ---------------------------------------------------
Y_unval             Column name of unvalidated outcome in the input
                    dataset.

Y                   Column name of validated outcome in the input
                    dataset. NAs in the input will be counted as
                    individuals not selected in phase two.

X_unval             Column names of unvalidated covariates in the input
                    dataset.

X                   Column names of validated covariates in the input
                    dataset. NAs in the input will be counted as
                    individuals not selected in phase two.

Z                   Column names of error-free covariates in the input
                    dataset.

Bspline             Column names of the B-spline basis in the input
                    dataset.

data                Dataset

hn_scale            Scale of the perturbation constant in the variance
                    estimation via the method of profile likelihood.
                    The default is 1.

noSE                Standard errors of the parameter estimators will
                    not be estimated when set to TRUE. The default is
                    FALSE.

TOL Convergence     The default is 0.0001.
criterion.          

MAX_ITER            Maximum number of iterations in the EM algorithm.
                    The default is 1000.

verbose             Print analysis details when set to TRUE. The
                    default is FALSE.
-----------------------------------------------------------------------

Table 2: Data dictionary for mock.vccc

-------------------------------------------------------------------------
Name           Status         Type         Description
-------------- -------------- ------------ ------------------------------
ID             error-free                  Patient ID

VL_unval       error-prone    continuous   Viral load (VL) at
                                           antiretroviral therapy (ART)
                                           initiation

VL_val         validated                   

ADE_unval      error-prone    binary       Had an AIDS-defining event
                                           (ADE) within one year of ART
                                           initiation: 1 - yes, 0 -- no

ADE_val        validated                   

CD4_unval      error-prone    continuous   CD4 count at ART initiation

CD4_val        validated                   

prior_ART      error-free     binary       Experienced ART before
                                           enrollment: 1 - yes, 0 - no

Sex            error-free     binary       Sex at birth of patient: 1 -
                                           male, 0 - female

Age            error-free     continuous   Age of patient
-------------------------------------------------------------------------
