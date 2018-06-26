\name{test.balsos}
\alias{test.balsos}
\title{Testing Measurement Level Assumptions.}
\description{Uses the algorithm discussed in Armstrong and Jacoby 2018) to test the intervality of a variable scaled by the Bayesian ALSOS method. }
\usage{
test.balsos(orig, chains, N=NULL)
}
\arguments{
  \item{orig}{Original dependent variable. }
  \item{chains}{N (obs) X K (iterations) matrix of optimally scaled points, 
  \item{N}{Number of draws from sampler to use.  If \code{NULL}, all draws are used. }
}
\value{
    Printed output giving the Bayesian p-value evaluating the null that ther eis no interesting difference between the original values and the optimally scaled values. 
}
\author{Dave Armstrong}
