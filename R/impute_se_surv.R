#' Impute Missing precision information for meta analysis of survival rates
#'
#' \code{impute_se_surv} imputes missing precision information for meta analysis
#'  of survival rates based on the imputation method proposed by Maruo et al.
#'  (submitting) and calculates survival rates and their SEs on the transformed
#'  scale.
#'
#' @param data a data frame that may include \code{St}, \code{LCL},
#'   \code{UCL}, \code{n}, \code{nt}, \code{ne}, \code{p}, and \code{methodvar}.
#'   If \code{LCL} and \code{UCL} are missing, at least one of \code{nt},
#'   \code{ne}, or \code{p} must be provided. Among them, \code{nt} should be
#'   provided.
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
#'   This variable is a character variable which can take the values  "plain",
#'   "log", or "log-log". If there are studies for which CIs are not missing and
#'   the estimation method is known, it may be specified. Default is
#'   \code{NULL}.If not identified, one of the three estimation methods above is
#'   automatically selected based on the symmetry of the CI of the transformed
#'   measure.
#' @param symmety_lim a numeric value for the threshold for the symmetry of the
#'   confidence interval, defined as the absolute value of the ratio of the
#'   distances from both confidence limits to the point estimate minus 1,
#'   on the selected transformed scale. If the symmetry exceeds this threshold,
#'   the confidence interval estimation method is not considered to be “plain,”
#'   “log,” or “log-log” and is replaced with \code{NA}.
#'   The default is 0.1.
#' @param percent_scale a logical value (TRUE or FALSE).
#'   If TRUE, St, LCL, UCL, and p are provided as percentages (0–100%).
#'   If FALSE, they are provided as proportions (0–1). The default is FALSE.
#'
#' @details In meta-analyses of survival rates, precision information (i.e.,
#' confidence interval) may be missing in some studies. Excluding studies
#' with missing precision information may reduce the precision and accuracy
#' of survival rate estimators in meta-analyses significantly.
#' This function imputes precision information using information commonly
#' available from study literature, such as sample size, number of events,
#' and risk set size at a time point of interest.
#' It then returns a data frame with the variables necessary for meta-analysis
#' added. Even when no data are missing, it can be used as a
#' preprocessing step in the meta-analysis of survival rates.
#'
#' @return a data frame with the following variables added to the input
#' data frame (\code{data}). These variables are used in a meta analysis.
#' Transformed scale is chosen with a \code{method} argument.
#' \describe{
#'   \item{\code{tr_St}}{survival rate on the transformed scale.}
#'   \item{\code{tr_SE}}{SE of survival rate on the transformed scale.}
#'   \item{\code{prt}}{parameter used for imputation. For detail, see Maruo
#'   et al. (submitting)}
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

impute_se_surv <- function(data, St, LCL, UCL, n, nt = NULL, ne = NULL,
                           p = NULL, method = "log-log", methodvar = NULL,
                           symmety_lim = 0.1, percent_scale = FALSE){
  if (!is.data.frame(data)){
    stop("data must be data.frame.")
  }
  data <- as.data.frame(data)
  notna <- !is.na(data[, St])
  nan <- nrow(data) - sum(notna)
  if (nan > 0) {
    warning(paste("St for", nan, "studies are missing. These are not imputed.")
    )
  }
  R0 <- sum(notna)
  Sta <- data[notna, St]
  LCLa <- data[notna, LCL]
  UCLa <- data[notna, UCL]
  na <- data[notna, n]
  if (sum(is.na(UCLa)) != 0 & is.null(nt) & is.null(p) & is.null(ne)){
    stop("If there are missing measurements in the LCL and UCL, then nt,
         ne or p must be specified.")
  }
  if (percent_scale) {
    Sta <- Sta / 100
    LCLa <- LCLa / 100
    UCLa <- UCLa / 100
  }
  if (is.null(nt)) {
    nta <- numeric(R0) * NA
  } else {
    nta <- data[notna, nt]
  }
  if (!is.null(ne)) {
    nea <- data[notna, ne]
    if (is.null(p)) {
      pa <- nea / na
    }
  }
  if (is.null(ne)){
    nea <- numeric(R0) * NA
    if (is.null(p)) {
      pa <- nea
    }
  }
  if (!is.null(p)) {
    pa <- data[notna, p]
    if (percent_scale) {
      pa <- pa / 100
    }
  }
  est <- c()
  Sta0n <- sum(Sta == 0, na.rm = TRUE)
  Sta1n <- sum(Sta == 1, na.rm = TRUE)
  if (Sta1n > 0) {
    warning(paste("St for", Sta1n,
                  "studies are 1 (100%), so they are replaced by 0.99.
                  Performance with respect to this replacement is not ensured.")
    )
  }
  if (Sta0n > 0) {
    warning(paste("St for", Sta0n,
                  "studies are 0, so they are replaced by 0.01.
                  Performance with respect to this replacement is not ensured.")
    )
  }
  Sta[Sta == 1] <- 0.99
  UCLa[Sta == 1] <- NA
  LCLa[Sta == 1] <- NA
  if (sum(Sta > 1 & Sta < 0) > 0) {
    stop("St (or St/100) must be in [0, 1].")
  }
  if (sum(LCLa > 1 & LCLa < 0, na.rm = TRUE) > 0) {
    stop("LCL (or LCL/100) must be in [0, 1].")
  }
  if (sum(UCLa > 1 & UCLa < 0, na.rm = TRUE) > 0) {
    stop("UCL (or UCL/100) must be in [0, 1].")
  }
  if (sum(pa > 1 & pa < 0, na.rm = TRUE) > 0) {
    stop("p (or p/100) must be in [0, 1].")
  }
  if (sum(nea > na, na.rm = TRUE) > 0) {
    stop("ne must be less than or equal to n")
  }
  if (sum(nta > na, na.rm = TRUE) > 0) {
    stop("nt must be less than or equal to n")
  }

  notna2 <- !(is.na(LCLa) & is.na(nta) & is.na(pa))
  R <- sum(notna2)
  Sta <- Sta[notna2]
  LCLa <- LCLa[notna2]
  UCLa <- UCLa[notna2]
  na <- na[notna2]
  nea <- nea[notna2]
  nta <- nta[notna2]
  pa <- pa[notna2]

  if (R0 - R > 0) {
    warning(paste("All of UCL, LCL, nt, ne, and p are missing in the",
                  R0 - R, "studies, so these are not imputed.")
    )
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
  if (method == "log-log") {
    Sta_t <- log(-log(Sta))
  }
  if (method == "log") {
    Sta_t <- log(Sta)
  }
  if (method == "plain"){
    Sta_t <- Sta
  }
  prta <- numeric(R)
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
        if (min(method_j1) > symmety_lim) {
          methodr <- NA
        }
      }
      if (!is.na(methodr)){
        if (methodr == "log-log") {
          eta_sqrt <- (log(-log(LCLr)) - log(-log(UCLr)))/(2 * qnorm(0.975)) *
            abs(log(Srt))
        }
        if (methodr == "log") {
          eta_sqrt <- (log(UCLr) - log(LCLr))/(2 * qnorm(0.975))
        }
        if (methodr == "plain") {
          eta_sqrt <- (UCLr - LCLr)/(2 * qnorm(0.975) *
                                       Srt)
        }
        methoda[r] <- methodr
      }
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
        if (!is.na(nrt)){
          prt <- eta_sqrt ^ 2 /
            sum(1 / (nrt:(nr - 1)) * 1 / ((nrt - 1):(nr - 2)))
          prta[r] <- prt
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
      prta[r] <- prt

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
  data2$tr_St <- NA
  data2$tr_St[notna][notna2] <- est[, 1]
  data2$tr_SE <- NA
  data2$tr_SE[notna][notna2] <-est[, 2]
  data2$prt <- NA
  data2$prt[notna][notna2] <- prta
  data2$imputed <- NA
  data2$imputed[notna][notna2] <- imputed
  return(data2)
}
