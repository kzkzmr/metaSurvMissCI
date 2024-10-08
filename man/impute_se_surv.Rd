% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_se_surv.R
\name{impute_se_surv}
\alias{impute_se_surv}
\title{Impute Missing precision information for meta analysis of survival rates}
\usage{
impute_se_surv(
  data,
  St,
  LCL,
  UCL,
  n,
  nt,
  ne,
  p,
  method = "log-log",
  methodvar = NULL
)
}
\arguments{
\item{data}{a data frame that may include \code{St}, \code{LCL},
\code{UCL}, \code{n}, \code{nt}, \code{ne}, \code{p}, and \code{methodvar}.
If \code{LCL} and \code{UCL} are missing, either \code{nt}, \code{ne}, or
\code{p} must be non-missing.}

\item{St}{a character value for name of t year survival rate variable
included in \code{data}.}

\item{LCL}{a character value for name of lower confidence limit for
survival rate variable included in \code{data}.}

\item{UCL}{a character value for name of upper confidence limit for
survival rate variable included in \code{data}.}

\item{n}{a character value for name of sample size variable included in
\code{data}.}

\item{nt}{a character value for name of variable for risk set size at
time point t included in \code{data}.}

\item{ne}{a character value for name of number of events variable included in
\code{data}.}

\item{p}{a character value for name of proportion of event occurrence
included in \code{data}.}

\item{method}{a character value specifying the estimation method for
confidence interval. Choose from "plain", "log", and "log-log".
Default is "log-log".}

\item{methodvar}{a character value for name of variable specifying the
estimation method for confidence interval for each row in \code{data}.
This variable is a character variable which can be "plain", "log", or
"log-log". If there are studies for which CIs are not missing and
the estimation method is known, it may be specified. Default is
\code{NULL}.If not identified, one of the three estimation methods above is
automatically selected based on the symmetry of the CI of the transformed
measure.}
}
\value{
a data frame with the following variables added to the input
data frame (\code{data}). These variables are used in a meta analysis.
Transformed scale is chosen with a \code{method} argument.
\describe{
  \item{\code{tr_St}}{survival rate on the transformed scale.}
  \item{\code{tr_SE}}{SE of survival rate on the transformed scale.}
  \item{\code{imputed}}{1: imputed, 0: not imputed.}
}
}
\description{
\code{impute_se_surv} imputes missing precision information for meta analysis
 of survival rates based on the imputation method proposed by Maruo et al.
 (submitting) and calculate survival rates and their SEs on the transformed
 scale.
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
data("metadata_chordoma")
mets <- impute_se_surv(data = metadata_chordoma, St = "PFS5y",
                       LCL = "PFSL5y", UCL = "PFSU5y", n = "n",
                       nt = "n_5yPFS", ne = "ne_PFS", p = "pr_PFS")

}
\references{
\itemize{
  \item Maruo K, Yamaguchi Y, Ishii R, Gosho M (submitting). Meta-analysis of
   survival rates when precision information is missing.
}
}
\seealso{
\code{\link{meta}}
}
