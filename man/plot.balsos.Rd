\name{plot.balsos}
\alias{plot.balsos}
\title{Plot Results from BALSOS}
\description{Plots the optimally scaled points with posterior 95\% credible intervals.}
\usage{
  \method{plot}{balsos}(x, ..., freq=TRUE, offset=.1)
}
\arguments{
  \item{x}{Object of class \code{balsos}.  }
  \item{freq}{Logical indicating whether you want the frequentist result plotted alongside the Bayesian result.}
  \item{offset}{If \code{freq=T}, the Bayesian points will be plotted at \code{x-offset} and the frequentist points will be plotted at \code{x+offset}.} 
  \item{...}{Other arguments to be passed down, currently not implement and may conflict with the lattice figure.  To change the figure the best advice would be to save the plot as an oject and use the \code{update} function to change its features. }
}
\value{
    A lattice graph produce by a call to \code{xyplot}.
}
\author{Dave Armstrong}
