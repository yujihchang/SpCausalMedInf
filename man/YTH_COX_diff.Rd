\name{YTH_COX_diff}
\alias{YTH_COX_diff}
\title{Semiparametric separable direct and indirect effects via Cox components}
\description{
Computes separable direct (DE) and indirect (IE) effects for semicompeting risks
with a time-to-event mediator (\code{X1}) and terminal outcome (\code{X2}),
allowing left truncation.
}
\usage{
YTH_COX_diff(time, A, X1, X2, D, X,
             expor = c(za1, za2, zb1, zb2),
             center = "TRUE", setting = NULL,
             EFFECT = c("DE","IE"), saving = "FALSE")
}
\arguments{
  \item{time}{Numeric vector of time points at which effects are evaluated.}
  \item{A}{Left-truncation times (length \eqn{n}); use 0 if none.}
  \item{X1}{Mediator event times (length \eqn{n}).}
  \item{X2}{Outcome event times (length \eqn{n}). Must satisfy \eqn{X1 \le X2}.}
  \item{D}{Event indicator for \code{X2} (1=event, 0=censor).}
  \item{X}{Covariate matrix of size \eqn{n \times p}.}
  \item{expor}{Exposure levels \code{c(za1, za2, zb1, zb2)}; e.g., \code{c(0,1,0,1)} or \code{c(1,2,1,2)}.}
  \item{center}{\code{"TRUE"} uses means; \code{"SET"} uses \code{setting}.}
  \item{setting}{Numeric vector of length \eqn{p-1} used when \code{center="SET"}.}
  \item{EFFECT}{Character vector: \code{"DE"}, \code{"IE"}, or both.}
  \item{saving}{Internal flag (\code{"TRUE"} or \code{"FALSE"}).}
}
\value{
A list with elements: \code{DE}, \code{IE}, \code{DE.sd}, \code{IE.sd},
and (when both effects are requested) cumulative hazard components
\code{LA_Sz0_za2}, \code{LA_Sz1_za2}, \code{LA_Sz0_za1}, \code{LA_Sz1_za1}.
}
\examples{
## See README for a full example.
}
