#' Mock VCCC dataset.
#'
#' A simulated dataset constructed to imitate the Vanderbilt Comprehensive Care
#' Clinic (VCCC) patient records, which have been fully validated and therefore
#' contain validated and unvalidated versions of all variables. The VCCC
#' cohort is a good candidate for the purpose of illustration. The
#' data presented in this section are a mocked-up version of the actual data
#' due to confidentiality, but the data structure and features, such as mean
#' and variability, closely resemble the real dataset.
#'
#' @format A data frame with 2087 rows and 8 variables:
#' \describe{
#'   \item{ID}{patient ID}
#'   \item{VL_unval}{viral load at antiretroviral therapy (ART) initiation,
#'   error-prone outcome, continuous}
#'   \item{VL_val}{viral load at antiretroviral therapy (ART) initiation,
#'   validated outcome, continuous}
#'   \item{ADE_unval}{having an AIDS-defining event (ADE) within one year
#'   of ART initiation, error-prone outcome, binary}
#'   \item{ADE_val}{having an AIDS-defining event (ADE) within one year of ART
#'   initiation, validated outcome, binary}
#'   \item{CD4_unval}{CD4 count at ART initiation, error-prone covariate,
#'   continuous}
#'   \item{CD4_val}{CD4 count at ART initiation, validated covariate,
#'   continuous}
#'   \item{prior_ART}{whether patient is ART naive at enrollment, error-free
#'   covariate, binary}
#'   \item{Sex}{sex of patient, 1 indicates male and 0 indicates
#'   female & error-free covariate, binary}
#'   \item{Age}{age of patient, error-free covariate, continuous}
#' }
#' @source \url{http://www.diamondse.info/}
"mock.vccc"
