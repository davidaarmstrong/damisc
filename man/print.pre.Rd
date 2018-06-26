\name{print.pre}
\alias{print.pre}
\title{Print method for objects of class pre}
\description{Prints the output from an object of class \code{pre}.  The function prints all components of the calculation and optionally simulated confidence bounds. }
\usage{
\method{print}{pre}(x, ..., sim.ci=.95)
}
\arguments{
  \item{x}{An object of class \code{pre}.}
  \item{sim.ci}{Coverage for the simulated confidence interval, if \code{sim=TRUE} in the call to \code{pre}.}
  \item{...}{Other arguments passed to print, currently not implemented}
}
\author{Dave Armstrong
}
\seealso{\code{pre}
}
