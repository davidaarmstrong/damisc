\name{secondDiff}
\alias{secondDiff}
\title{Calculate Cross-Derivative and its Variability}
\description{Calculates the cross-derivative required to evaluate interactions in logistic/probit regression models.}
\usage{
secondDiff(obj, vars, data, method=c("AME", "MER"), typical = NULL)
}
\arguments{
  \item{obj}{An object of class \code{glm} that will be used to find the cross-derivative. }
  \item{vars}{A vector of two variables to be used in calculating the derivative. }
  \item{data}{A data frame.}
  \item{method}{Indicate whether you want to use average marginal effects (AME) or marginal effects at representative values (MER). }
  \item{typical}{A named vector of values at which to hold variables constant. }
}
\value{
A data frame with standardized quantitative variables
}
\author{Dave Armstrong}
