\name{test.balsos}
\alias{test.balsos}
\title{Testing Measurement Level Assumptions.}
\description{Uses the algorithm discussed in Armstrong and Jacoby 2018) to test the intervality of a variable scaled by the Bayesian ALSOS method. }
\usage{
test.balsos(obj, N=NULL)
}
\arguments{
  \item{obj}{An object of class \code{balsos}.}
  \item{N}{Number of draws from sampler to use.  If \code{NULL}, all draws are used. }
}
\value{
    Printed output giving the Bayesian p-value evaluating the null that ther eis no interesting difference between the original values and the optimally scaled values. 
}
\author{Dave Armstrong}
