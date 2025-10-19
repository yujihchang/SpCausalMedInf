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
\dontrun{
## --- Example using the bundled dataset expdata ---------------------------
## If you placed expdata.rda under the package's data/ folder, you can:
data("expdata", package = "SpCausalMedInf")

## Use the dataset
work_ttem <- expdata

## Add a tiny jitter to avoid tied times
set.seed(123)
ee <- runif(nrow(work_ttem), 0, 0.001)

## Map columns (adjust names to your actual data)
X1 <- work_ttem$hcc.time + ee
X2 <- work_ttem$die.time + ee
D  <- work_ttem$dieall.case
X  <- as.matrix(work_ttem$HBSAG)  ## if multiple covariates: as.matrix(work_ttem[, c("HBSAG","...")])
A  <- work_ttem$AGE
t  <- sort(X2)

## Run
NPMLE <- YTH_COX_diff(t, A, X1, X2, D, X,
                      expor = c(0,1,0,1),
                      center = "TRUE", setting = NULL,
                      EFFECT = c("DE","IE"))
# center = "TRUE": set covariates to their sample means.

NPMLE_set <- YTH_COX_diff(t, A, X1, X2, D, X,
                          expor = c(0,1,0,1),
                          center = "SET", setting = c(40, 1, 0),
                          EFFECT = c("DE","IE"))
# center = "SET": fix covariates at specified levels.
# Here: age = 40, smoker = 1 (yes), alcohol = 0 (no).
# Note: 'setting' supplies values for X[, -1] in that column order.
}
}

