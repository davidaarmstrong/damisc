\name{testLoess}
\alias{testLoess}
\title{Significance Test for Loess vs. LM}
\description{Calculates an F test to evaluate significant differences between a LOESS model and a parametric alternative estimated with \code{lm}}
\usage{
testLoess(lmobj, loessobj, alpha=.05)
}
\arguments{
  \item{lmobj}{An object of class \code{lm}.}
  \item{loessobj}{An object of class \code{loess}.}
  \item{alpha}{Desired Type I error rate of test.}
}
\value{
    Printed output describing the results of the test.
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
linmod <- lm(prestige ~ income, data=Prestige)
lomod <- loess(prestige ~ income, data=Prestige)
testLoess(linmod, lomod)
}
