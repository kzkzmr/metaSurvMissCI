#' Forestplot for meta analysis of survival rates
#'
#' \code{forest_surv} is a function that draws a forest plot using a
#' \code{metagen} object of a survival rate meta-analysis based on the
#' \code{metagen} function. It is assumed that the data frame output from
#' the \code{impute_se_surv} function is used in the \code{metagen} function.
#' @param metares a \code{metagen} object for the meta analysis of survival
#'   rate.
#' @param method a character value specifying the estimation method for
#'   confidence interval used in the \code{impute_se_surv} function.
#'   Choose from "plain", "log", and "log-log". Default is "log-log".
#' @param percent a logical value that specifies the scale of the display;
#'   \code{TRUE} for a percentage scale (0 to 100),
#'   \code{FALSE} for a rate scale (0 to 1).
#'   Default is \code{TRUE}.
#' @param xlim a numeric vector of length 2 that specifies horizontal limits of
#'   the plot region. If unspecified, the function sets the horizontal plot
#'   limits to some sensible values.
#' @param estlab a character value for labels for estimation target.
#'   Default is "Survival (\%)".
#' @param digits a numeric value specifying the number of decimal places to
#'   be displayed. Default is \code{1}.
#' @param ... additional arguments to be passed to the \code{forest} function.
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
#' @return a forest plot created with the \code{forest} function in the
#'   \code{meta} package.
#'
#'
#' @seealso \code{\link{meta}}
#'
#' @examples
#' library(metaSurvMissCI)
#' library(meta)
#' data("metadata_chordomas")
#' impdata <- impute_se_surv(data = metadata_chordomas, St = "PFS5y",
#'                        LCL = "PFSL5y", UCL = "PFSU5y", n = "n",
#'                        nt = "n_5yPFS", ne = "ne_PFS", p = "pr_PFS")
#' meta_imp <- metagen(TE = tr_St, seTE = tr_SE, studlab = Study,
#'                     data = impdata, method.tau = "REML")
#' forest_surv(meta_imp, xlim = c(30, 100), estlab = "5-year PFS (%)")
#'
#' @importFrom meta forest
#'
#' @export

forest_surv <- function(metares, method = "log-log", percent = TRUE,
                          xlim = NULL, estlab = "Survival (%)", digits = 1,
                          ...){
  metar_f <- metares
  cc <- 1
  if (percent) {
    cc <- 100
  }
  if (method == "log-log"){
    metar_f$TE <- exp(-exp(metares$TE)) * cc
    metar_f$lower <- exp(-exp(metares$upper)) * cc
    metar_f$upper <- exp(-exp(metares$lower)) * cc
    metar_f$TE.random <- exp(-exp(metares$TE.random)) * cc
    metar_f$lower.random <- exp(-exp(metares$upper.random)) * cc
    metar_f$upper.random <- exp(-exp(metares$lower.random)) * cc
    metar_f$TE.common <- exp(-exp(metares$TE.common)) * cc
    metar_f$lower.common <- exp(-exp(metares$upper.common)) * cc
    metar_f$upper.common <- exp(-exp(metares$lower.common)) * cc
  }
  if (method == "log"){
    metar_f$TE <- exp(metares$TE) * cc
    metar_f$lower <- exp(metares$lower) * cc
    metar_f$upper <- exp(metares$upper) * cc
    metar_f$TE.random <- exp(metares$TE.random) * cc
    metar_f$lower.random <- exp(metares$lower.random) * cc
    metar_f$upper.random <- exp(metares$upper.random) * cc
    metar_f$TE.common <- exp(metares$TE.common) * cc
    metar_f$lower.common <- exp(metares$lower.common) * cc
    metar_f$upper.common <- exp(metares$upper.common) * cc
  }
  if (method == "plain"){
    metar_f$TE <- metares$TE * cc
    metar_f$lower <- metares$lower * cc
    metar_f$upper <- metares$upper * cc
    metar_f$TE.random <- metares$TE.random * cc
    metar_f$lower.random <- metares$lower.random * cc
    metar_f$upper.random <- metares$upper.random * cc
    metar_f$TE.common <- metares$TE.common * cc
    metar_f$lower.common <- metares$lower.common * cc
    metar_f$upper.common <- metares$upper.common * cc
  }
  forest(metar_f, xlim = xlim, leftcols = c("studlab", "n"),
         leftlabs = c("Study                   ", "N"),
         rightlabs = c(estlab, "[95%CI]", "W(common)", "W(random)"),
         digits = digits, ...)
}
