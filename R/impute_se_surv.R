#' Impute Missing precision information for meta analysis of survival rates
#'
#' \code{impute_se_surv} imputes missing precision information for meta analysis
#'  of survival rates based on the imputation method proposed by Maruo et al.
#'  (submitting) and calculate survival rates and their SEs on the transformed
#'  scale.
#'
#' @param data a data frame that may include \code{St}, \code{LCL},
#'   \code{UCL}, \code{n}, \code{nt}, \code{ne}, \code{p}, and \code{methodvar}.
#'   If \code{LCL} and \code{UCL} are missing, either \code{nt}, \code{ne}, or
#'   \code{p} must be non-missing.
#' @param St a character value for name of t year survival rate variable
#'   included in \code{data}.
#' @param LCL a character value for name of lower confidence limit for
#'   survival rate variable included in \code{data}.
#' @param UCL a character value for name of upper confidence limit for
#'   survival rate variable included in \code{data}.
#' @param n a character value for name of sample size variable included in
#'   \code{data}.
#' @param nt a character value for name of variable for risk set size at
#'   time point t included in \code{data}.
#' @param ne a character value for name of number of events variable included in
#'   \code{data}.
#' @param p a character value for name of proportion of event occurrence
#'   included in \code{data}.
#' @param method a character value specifying the estimation method for
#'   confidence interval. Choose from "plain", "log", and "log-log".
#'   Default is "log-log".
#' @param methodvar a character value for name of variable specifying the
#'   estimation method for confidence interval for each row in \code{data}.
#'   This variable is a character variable which can be "plain", "log", or
#'   "log-log". If there are studies for which CIs are not missing and
#'   the estimation method is known, it may be specified. Default is
#'   \code{NULL}.If not identified, one of the three estimation methods above is
#'   automatically selected based on the symmetry of the CI of the transformed
#'   measure.
#'
#' @details In meta-analyses of survival rates, precision information (i.e.,
#' confidence interval) may be missing in some studies. Excluding studies
#' with missing precision information may reduce the precision and accuracy
#' of survival rate estimators in meta-analyses significantly.
#' This function impute precision information using information commonly
#' available from study literature, such as sample size, number of events,
#' and risk set size at a time point of interest.
#' It then returns a data frame with the variables necessary for meta-analysis
#' added. Even in the absence of missing data, it can be used as a
#' pre-processing step in the meta-analysis of survival rates.
#'
#' @return a data frame with the following variables added to the input
#' data frame (\code{data}). These variables are used in a meta analysis.
#' Transformed scale is chosen with a \code{method} argument.
#' \describe{
#'   \item{\code{tr_St}}{survival rate on the transformed scale.}
#'   \item{\code{tr_SE}}{SE of survival rate on the transformed scale.}
#'   \item{\code{imputed}}{1: imputed, 0: not imputed.}
#' }
#'
#' @references \itemize{
#'   \item Maruo K, Yamaguchi Y, Ishii R, Gosho M (submitting). Meta-analysis of
#'    survival rates when precision information is missing.
#' }
#'
#' @seealso \code{\link{meta}}
#'
#' @examples
#' library(metaSurvMissCI)
#' data("metadata_chordoma")
#' mets <- impute_se_surv(data = metadata_chordoma, St = "PFS5y",
#'                        LCL = "PFSL5y", UCL = "PFSU5y", n = "n",
#'                        nt = "n_5yPFS", ne = "ne_PFS", p = "pr_PFS")
#'
#' @importFrom stats qnorm
#'
#' @export

impute_se_surv <- function(data, St, LCL, UCL, n, nt, ne, p,
                           method = "log-log", methodvar = NULL){
  if (!is.data.frame(data)){
    stop("data must be data.frame.")
  }
  R <- nrow(data)
  Sta <- data[, St]
  LCLa <- data[, LCL]
  UCLa <- data[, UCL]
  na <- data[, n]
  nta <- data[, nt]
  nea <- data[, ne]
  pa <- data[, p]
  est <- c()
  Sta[Sta == 0] <- 0.01
  UCLa[Sta == 0] <- NA
  LCLa[Sta == 0] <- NA
  Sta[Sta == 1] <- 0.99
  UCLa[Sta == 1] <- NA
  LCLa[Sta == 1] <- NA
  if (sum(Sta > 1 & Sta < 0) > 0) {
    stop("St must be in [0, 1].")
  }
  if (sum(LCLa > 1 & LCLa < 0, na.rm = TRUE) > 0) {
    stop("LCL must be in [0, 1].")
  }
  if (sum(UCLa > 1 & UCLa < 0, na.rm = TRUE) > 0) {
    stop("UCL must be in [0, 1].")
  }
  if (sum(pa > 1 & pa < 0, na.rm = TRUE) > 0) {
    stop("p must be in [0, 1].")
  }
  if (sum(nea > na, na.rm = TRUE) > 0) {
    stop("ne must be less than or equal to n")
  }
  if (sum(nta > na, na.rm = TRUE) > 0) {
    stop("nt must be less than or equal to n")
  }

  imputed <- numeric(R)
  if (is.null(methodvar)) {
    methoda <- numeric(R) * NA
  } else {
    methoda <- data[, methodvar]
  }
  method_lst <- c("log-log", "log", "plain")
  if (!(method %in% method_lst)) {
    stop('method must be selected from c("log-log", "log", "plain")')
  }
  if (sum(is.na(Sta)) > 0) {
    stop("There are missing values in the survival rate variable (St).")
  }
  if (method == "log-log") {
    Sta_t <- log(-log(Sta))
  }
  if (method == "log") {
    Sta_t <- log(Sta)
  }
  if (method == "plain"){
    Sta_t <- Sta
  }
  for (r in 1:R){
    Srt <- Sta[r]
    UCLr <- UCLa[r]
    LCLr <- LCLa[r]
    nr <- na[r]
    nre <- nea[r]
    nrt <- nta[r]
    pr <- pa[r]
    Srt_t <- Sta_t[r]
    methodr <- methoda[r]

    if (!is.na(LCLr) & !is.na(UCLr)) {
      if (!is.na(methodr)) {
        if (!(methodr %in% method_lst)) {
          warning('method must be selected from c("log-log", "log", "plain")
                  method is replaced by NA.')
          methodr <- NA
        }
      }
      if (is.na(methodr)) {
        method_j1 <- abs(c((log(-log(LCLr)) - log(-log(Srt))) /
                             (log(-log(Srt)) - log(-log(UCLr))),
                           (log(UCLr) - log(Srt)) / (log(Srt) - log(LCLr)),
                           (UCLr - Srt) / (Srt - LCLr)) - 1)
        methodr <- method_lst[which.min(method_j1)]
        if (min(method_j1) > 0.1) {
          methodr <- NA
        }
      }
      if (methodr == "log-log"){
        eta_sqrt <- (log(-log(LCLr)) - log(-log(UCLr))) / (2 * qnorm(0.975)) *
          abs(log(Srt))
      }
      if (methodr == "log"){
        eta_sqrt <- (log(UCLr) - log(LCLr)) / (2 * qnorm(0.975))
      }
      if (methodr == "plain"){
        eta_sqrt <- (UCLr - LCLr) / (2 * qnorm(0.975) * Srt)
      }
      methoda[r] <- methodr
      if (is.na(methodr)) {
        LCLr <- NA
        UCLr <- NA
        SEr <- NA
      } else {
        if (method == "log-log"){
          SEr <- eta_sqrt / abs(log(Srt))
        }
        if (method == "log"){
          SEr <- eta_sqrt
        }
        if (method == "plain") {
          SEr <- eta_sqrt * Srt
        }
      }
    }

    if ((is.na(LCLr) | is.na(UCLr))){
      imputed[r] <- 1
      if (is.na(nr)){
        stop("n must be given when LCL or UCL are not given.")
      }
      if (is.na(nrt)){
        if (is.na(nre) & is.na(pr)){
          stop("ne or p must be given when LCL, UCL, and nt are not given.")
        }
        if (!is.na(nre)) {
          pr <- nre / nr
        }
        nrt <- nr * ((pr + 1) * Srt + pr - 1)/((pr - 1) * Srt + pr + 1)
      }
      nrt <- max(c(2, nrt))
      nrt <- min(c(nrt, nr - 1))
      prt <- ((nr + nrt) * Srt - nrt - nr) / ((nrt - nr) * Srt + nrt - nr)
      prt[prt > 1] <- 1
      prt[prt < 0] <- 0.01

      eta <- sum(1 / (nrt:(nr - 1)) * 1 / ((nrt - 1):(nr - 2))) * prt
      if (method == "log-log"){
        SEr <- sqrt(eta) / abs(log(Srt))
      }
      if (method == "log"){
        SEr <- sqrt(eta)
      }
      if (method == "plain") {
        SEr <- sqrt(eta) * Srt
      }
    }
    est <- rbind(est, t(c(Srt_t, SEr)))
  }
  data2 <- data
  data2$tr_St <- est[, 1]
  data2$tr_SE <- est[, 2]
  data2$imputed <- imputed
  return(data2)
}
