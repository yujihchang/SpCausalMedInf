\name{YTH_COX_diff}
\alias{YTH_COX_diff}
\title{Nonparametric separable direct and indirect effects via hazard NPMLE}
\description{
Computes separable DE and IE effects using a nonparametric hazard-based approach
for semicompeting risks data with possible left truncation.
}
\usage{
YTH_haz_complete(t, A, X1, X2, D, Za, expo = c(za1,za2,zb1,zb2), EFFECT = c("DE","IE"))
}
\arguments{
  \item{t}{Numeric vector of time points at which effects are evaluated.}
  \item{A}{Left-truncation times (length \eqn{n}); use 0 if none.}
  \item{X1}{Event time for the mediative process (length \eqn{n}).}
  \item{X2}{Event time for the terminal outcome (length \eqn{n}). Must satisfy \eqn{X1 \le X2}.}
  \item{D}{Event indicator for \code{X2} (1=event, 0=censor).}
  \item{Za}{Numeric exposure vector (\eqn{n \times 1}).}
  \item{expo}{Exposure levels \code{c(za1, za2, zb1, zb2)}; e.g., \code{c(0,1,0,1)} or \code{c(1,2,1,2)}.}
  \item{EFFECT}{Character vector: \code{"DE"}, \code{"IE"}, or both.}
}
\value{
A list with elements: \code{DE}, \code{IE}, \code{DE.sd}, \code{IE.sd}.
}
\seealso{\code{\link{YTH_COX_diff}}}
\examples{
## See README for a full example.
}

