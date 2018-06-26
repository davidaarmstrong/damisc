\name{panel.ci}
\alias{panel.ci}
\title{Lattice panel function for confidence intervals}
\description{This panel function is defined to plot confidence intervals in a multi-panel lattice display. Note, both lower and upper must be passed directly to \code{xyplot} as they will be passed down to the prepanel function.}
\usage{
	panel.ci(x,y,subscripts, lower,upper,zl)
}
\arguments{
  \item{x,y}{Data from the call to \code{xyplot}.}
  \item{subscripts}{Variable used to created the juxtaposed panels.}
  \item{lower, upper}{95\% lower and upper bounds of \code{y}.}
  \item{zl}{Logical indicating whether or not a horizontal dotted line at zero is desired.}
}
\author{Dave Armstrong}
\examples{
library(car)
library(lattice)
data(Ornstein)
mod <- lm(interlocks ~ log(assets)*nation, data=Ornstein)
mod.out <- intQualQuant(mod, c("log(assets)", "nation"), 
	n=25, plot=FALSE, type="facs")
xyplot(fit ~ x | contrast, data=mod.out, 
	xlab = "Assets", ylab = "Difference In Fitted Values", 
	lower=mod.out$lower, upper=mod.out$upper, zl=TRUE, 
	prepanel=prepanel.ci, panel=panel.ci)
}

