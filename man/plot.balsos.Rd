\name{plot.balsos}
\alias{plot.balsos}
\title{Plot Results from BALSOS}
\description{Plots the optimally scaled points with posterior 95\% credible intervals.}
\usage{
plot.balsos(x, ..., freq=TRUE, offset=.1)
}
\arguments{
  \item{x}{Object of class \code{balsos}.  }
  \item{freq}{Logical indicating whether you want the frequentist result plotted alongside the Bayesian result.}
  \item{offset}{If \code{freq=T}, the Bayesian points will be plotted at \code{x-offset} and the frequentist points will be plotted at \code{x+offset}.} 
}
\value{
    A lattice graph produce by a call to \code{xyplot}.
}
\author{Dave Armstrong}
