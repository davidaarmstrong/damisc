\name{NKnots}
\alias{NKnots}
\title{AIC and BIC selection of number of spline knots}
\description{Calculates AIC and BIC for the selection of knots in a spline over values (potentially including polynomials) up to a user-defined maximum.}
\usage{
    NKnots(form, var, data, degree=3, min.knots=1, 
       max.knots=10, includePoly = FALSE, plot=FALSE, 
       criterion=c("AIC", "BIC"), K=10)
}
\arguments{
  \item{form}{A formula detailing the model for which smoothing is to be evaluated.}
  \item{var}{A character string identifying the variable for which smoothing is to be evaluated.}
  \item{data}{Data frame providing values of all variables in \code{form}.}
  \item{degree}{Degree of polynomial in B-spline basis functions.}
  \item{min.knots}{Minimum number of internal B-spline knots to be evaluated.}
  \item{max.knots}{Maximum number of internal B-spline knots to be evaluated.}
  \item{includePoly}{Include linear and polynomial models up to, and including \code{degree}-th order polynomials.}
  \item{plot}{Logical indicating whether a plot should be returned.}
  \item{criterion}{Statistical criterion to minimize in order to find the best number of knots.}
  \item{K}{If \code{criterion} is \sQuote{CV}, then \code{K} is the number of sub-samples used in the cross-validation}
}
\value{
    A plot, if \code{plot=TRUE}, otherwise a data frame with the degrees of freedom and corresponding fit measure. 
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
NKnots(prestige ~ education + type, var="income", data=na.omit(Prestige), plot=FALSE)
}
