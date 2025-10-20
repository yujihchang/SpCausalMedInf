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

  \item{X2}{Outcome event times (length \eqn{n}). Must satisfy \eqn{X1 \le X2}. 
  If \eqn{X2 > X1} (strictly), the mediator event is considered to have occurred. 
  If \eqn{X2 = X1}, the mediator is treated as \emph{censored} at \eqn{X1} (i.e., not observed). 
  If your data indicate the mediator occurred yet \eqn{X1 = X2}, add a very small noise to \eqn{X2} to separate them (e.g., replace \eqn{X2} by \eqn{X2 + \epsilon} with a tiny \eqn{\epsilon>0}).}

  \item{D}{Event indicator for \code{X2} (1 = event, 0 = censor).}

  \item{X}{An \eqn{n \times p} matrix. The \textbf{first column} must be the exposure variable. 
  The remaining \eqn{p-1} columns (if any) are covariates and must appear \textbf{after} the exposure in this order. 
  At minimum, \eqn{X} must include the exposure column; covariate columns are optional (i.e., \eqn{p=1} is allowed).}

  \item{expor}{Exposure levels \code{c(za1, za2, zb1, zb2)}; e.g., \code{c(0,1,0,1)} or \code{c(1,2,1,2)}.}

  \item{center}{Character flag: \code{"TRUE"} uses sample means of covariates; \code{"SET"} uses \code{setting}.}

  \item{setting}{Numeric vector of length \eqn{p-1} used when \code{center = "SET"} (values correspond to \eqn{X[, -1]} in column order).}

  \item{EFFECT}{Character vector: \code{"DE"}, \code{"IE"}, or both.}

  \item{saving}{Internal flag (\code{"TRUE"} or \code{"FALSE"}; default \code{"FALSE"}). 
  When \code{"FALSE"}, temporary objects created during each run are \emph{not} automatically cleared. 
  When \code{"TRUE"}, those temporaries are cleared after each run—this costs extra time but reduces peak memory usage, which can help with larger datasets.}
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
work_ttem <- Data

## Add a tiny jitter to avoid tied times
set.seed(123)
ee <- runif(nrow(work_ttem), 0, 0.001)

## Map columns (adjust names to your actual data)
X1 <- work_ttem$X1 + ee
X2 <- work_ttem$X2 + ee
D  <- work_ttem$D2
X  <- as.matrix(work_ttem[,-(1:4)])  ## if multiple covariates: as.matrix(work_ttem[, c("HBSAG","...")])
A  <- work_ttem$X2*0 ## if there is no truncation
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






