\name{ordfit}
\alias{ordfit}
\alias{print.ordfit}
\title{Fit Statistics for Proportional Odds Logistic Regression Models}
\description{For objects of class \code{polr}, it calculates a number of fit statistics and specification tests.}
\usage{
   ordfit(obj)
}
\arguments{
  \item{obj}{A model object of class \code{polr}.}
}
\value{
An object of class \code{ordfit} which is a matrix containing statistics and specification tests.
}
\references{
Lipsitz, S. R., Fitzmaurice, G. M. and Mohlenberghs, G.  1996.  Goodness-of-fit Tests for Ordinal Response Regression Models.  Applied Statistics, 45: 175-190.\cr
Pulkstenis, E. and Robinson, T. J. 2004.  Goodness-of-fit Test for Ordinal Response Regression Models.  Statistics in Medicine, 23: 999-1014. \cr
Fagerland, M. W. and Hosmer, D. W.  2013.  A Goodness-of-fit Test for the Proportional Odds Regression Model.  Statistics in Medicine 32(13): 2235-2249.
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(MASS)
data(france)
polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
ordfit(polr.mod)
}