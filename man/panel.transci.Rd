\name{panel.transci}
\alias{panel.transci}
\title{Lattice panel function for translucent confidence intervals}
\description{This panel function is defined to plot translucent confidence intervals in a single-panel, grouped (i.e., superposed) lattice display. Note, both lower and upper must be passed directly to \code{xyplot} as they will be passed down to the panel function.}
\usage{
	panel.transci(x,y,groups, lower,upper,...)
}
\arguments{
  \item{x,y}{Data from the call to \code{xyplot}.}
  \item{groups}{Variable used to created the superposed panels.}
  \item{lower, upper}{95\% lower and upper bounds of \code{y}.}
  \item{...}{Other arguments to be passed down to the plotting functions.}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
