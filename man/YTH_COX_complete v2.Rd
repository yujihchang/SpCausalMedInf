\name{YTH_haz_complete}
\alias{YTH_haz_complete}
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
  \item{A}{Left truncation times (length \(n\)); use 0 if none.}
  \item{X1}{Event time for the mediative process (length \(n\)).}
  \item{X2}{Event time for the terminal outcome (length \(n\)). \strong{Must satisfy} \(X1 \le X2\).}
  \item{D}{Event indicator for \code{X2} (1=event, 0=censor).}
  \item{Za}{Numeric exposure vector (\(n \times 1\)).}
  \item{expo}{Exposure levels \code{c(za1, za2, zb1, zb2)}; e.g., \code{c(0,1,0,1)} or \code{c(1,2,1,2)}.}
  \item{EFFECT}{Character vector choosing effects to compute: \code{"DE"}, \code{"IE"}, or both.}
}
\value{
A list with elements:
\item{DE}{Cumulative direct effect at \code{t}.}
\item{IE}{Cumulative indirect effect at \code{t}.}
\item{DE.sd, IE.sd}{Associated variability measures (variance scale).}
}
\seealso{\code{\link{YTH_COX_diff}}}
\examples{
## See README for a full example.
}
