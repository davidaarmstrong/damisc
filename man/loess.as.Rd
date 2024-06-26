% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{loess.as}
\alias{loess.as}
\alias{as.crit}
\alias{opt.span}
\title{Fit a local polynomial regression with automatic smoothing parameter selection}
\usage{
loess.as(
  x,
  y,
  degree = 1,
  criterion = c("aicc", "gcv"),
  family = c("gaussian", "symmetric"),
  user.span = NULL,
  plot = FALSE,
  ...
)

as.crit(x)

opt.span(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95))
}
\arguments{
\item{x}{An object of class \code{loess}.}

\item{y}{a vector of response values}

\item{degree}{the degree of the local polynomials to be used. It can ben 0, 1 or 2.}

\item{criterion}{The criterion used to find the optimal span}

\item{family}{if "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function.}

\item{user.span}{the user-defined parameter which controls the degree of smoothing.}

\item{plot}{if \code{TRUE}, the fitted curve or surface will be generated.}

\item{...}{control parameters.

Fit a local polynomial regression with automatic smoothing parameter selection. The predictor x can either one-dimensional or two-dimensional. This function was taken directly from `fANCOVA` version 0.5-1 and is wholly attributed to its author Xiao-Feng Wang.}

\item{model}{An object of class \code{loess}.}

\item{span.range}{The range in which to look for the optimal span}
}
\value{
An object of class \code{loess} with optimal smoothing parameter.
}
\description{
Fit a local polynomial regression with automatic smoothing parameter selection. Two methods are available for the selection of the smoothing parameter: bias-corrected Akaike information criterion (aicc); and generalized cross-validation (gcv).
}
\examples{
## Fit Local Polynomial Regression with Automatic Smoothing Parameter Selection
n1 <- 100
x1 <- runif(n1,min=0, max=3)
sd1 <- 0.2
e1 <- rnorm(n1,sd=sd1)
y1 <- sin(2*x1) + e1
(y1.fit <- loess.as(x1, y1, plot=TRUE))

}
\author{
X.F. Wang
}
