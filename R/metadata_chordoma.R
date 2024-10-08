#' Data for Systematic Review of Particle Beam Therapy for Chordoma
#'
#' The data from the systematic review of particle beam therapy for chordoma
#' (Saito et al., 2024). Data includes 7 studies for particle beam therapy
#' and 5-year progression-free survival (PFS).
#'
#' @name metadata_chordoma
#'
#' @docType data
#'
#' @format A data frame with 7 observations (studies) on the following 8
#' variables:
#'   \describe{
#'     \item{\code{Study}}{a character variable for study identifier.}
#'     \item{\code{PFS5y}}{a numeric variable for 5 year PFS.}
#'     \item{\code{PFSL5y}}{a numeric variable for the 95\% lower confidence
#'                          limit of 5 year PFS.}
#'     \item{\code{PFSU5y}}{a numeric variable for the 95\% upper confidence
#'                          limit of 5 year PFS.}
#'     \item{\code{n}}{an integer variable for the sample size.}
#'     \item{\code{ne_PFS}}{an integer variable for the number of PFS events.}
#'     \item{\code{n_5yPFS}}{an integer variable for the risk set size for PFS
#'                           at 5 year. All values are \code{NA}.}
#'     \item{\code{pr_PFS}}{a numeric variable for the proportion of event
#'                          occurrence for PFS. Since \code{ne_PFS} is missing
#'                          in two studies, they are imputed by the mean of
#'                          the other studies.}
#'   }
#'
#' @examples data(metadata_chordoma)
#'
#' @references \itemize{
#'   \item Saito T et al. Systematic review and meta-analysis of particle beam
#'   therapy versus photon radiotherapy for skull base chordoma: TRP-chordoma
#'   2024. \emph{Cancers}. 2024;16:2569.\doi{10.3390/cancers16142569}.
#'   }
NULL
