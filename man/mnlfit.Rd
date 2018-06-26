\name{mnlfit}
\alias{mnlfit}
\alias{print.mnlfit}
\title{Fit Statistics and Specification Test for Multinomial Logistic Regression}
\description{
Provides fit statistics (pseudo R-squared values) and the Fagerland, Hosmer and Bonfi (2008) specification test for Multinomial Logistic Regression models.}
\usage{
mnlfit(obj, permute = FALSE)
}
\arguments{
  \item{obj}{An object of class \code{multinom}}
  \item{permute}{Logical indicating whether to check all base categories for the Fagerland et. al. specification test.}
}
\value{
   A list with elements: 
   \item{result}{Fit statistics.}
   \item{permres}{The results of the base category permutation exercise.}
}
\references{
Fagerland M. W., D. W. Hosmer and A. M. Bonfi.  2008.  \sQuote{Multinomial goodness-of-fit tests for logistic regression models.}  Statistics in Medicine.  27(21): 4238-4253.
}
\author{Dave Armstrong}
\examples{
library(nnet)
data(france)
mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
mnlfit(mnl.mod)
}
