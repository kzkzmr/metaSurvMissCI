% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forest_surv.R
\name{forest_surv}
\alias{forest_surv}
\title{Forestplot for meta analysis of survival rates}
\usage{
forest_surv(
  metares,
  method = "log-log",
  percent = TRUE,
  xlim = NULL,
  estlab = "Survival (\%)",
  digits = 1,
  ...
)
}
\arguments{
\item{metares}{a \code{metagen} object for the meta analysis of survival
rate.}

\item{method}{a character value specifying the estimation method for
confidence interval used in the \code{impute_se_surv} function.
Choose from "plain", "log", and "log-log". Default is "log-log".}

\item{percent}{a logical value that specifies the scale of the display;
\code{TRUE} for a percentage scale (0 to 100),
\code{FALSE} for a rate scale (0 to 1).
Default is \code{TRUE}.}

\item{xlim}{a numeric vector of length 2 that specifies horizontal limits of
the plot region. If unspecified, the function sets the horizontal plot
limits to some sensible values.}

\item{estlab}{a character value for labels for estimation target.
Default is "Survival (\%)".}

\item{digits}{a numeric value specifying the number of decimal places to
be displayed. Default is \code{1}.}

\item{...}{additional arguments to be passed to the \code{forest} function.}
}
\value{
a forest plot created with the \code{forest} function in the
  \code{meta} package.
}
\description{
\code{forest_surv} is a function that draws a forest plot using a
\code{metagen} object of a survival rate meta-analysis based on the
\code{metagen} function. It is assumed that the data frame output from
the \code{impute_se_surv} function is used in the \code{metagen} function.
}
\details{
In meta-analyses of survival rates, precision information (i.e.,
confidence interval) may be missing in some studies. Excluding studies
with missing precision information may reduce the precision and accuracy
of survival rate estimators in meta-analyses significantly.
This function impute precision information using information commonly
available from study literature, such as sample size, number of events,
and risk set size at a time point of interest.
It then returns a data frame with the variables necessary for meta-analysis
added. Even in the absence of missing data, it can be used as a
pre-processing step in the meta-analysis of survival rates.
}
\examples{
library(metaSurvMissCI)
library(meta)
data("metadata_chordoma")
impdata <- impute_se_surv(data = metadata_chordoma, St = "PFS5y",
                       LCL = "PFSL5y", UCL = "PFSU5y", n = "n",
                       nt = "n_5yPFS", ne = "ne_PFS", p = "pr_PFS")
meta_imp <- metagen(TE = tr_St, seTE = tr_SE, studlab = Study,
                    data = impdata, method.tau = "REML")
forest_surv(meta_imp, xlim = c(30, 100), estlab = "5-year PFS (\%)")

}
\seealso{
\code{\link{meta}}
}
