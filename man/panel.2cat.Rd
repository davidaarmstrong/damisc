\name{panel.2cat}
\alias{panel.2cat}
\title{Lattice panel function for confidence intervals with capped bars}
\description{This panel function is defined to plot confidence intervals in a multi-panel lattice display where the x-variable is categorical. Note, both lower and upper must be passed directly to \code{xyplot} as they will be passed down to the panel function.}
\usage{
	panel.2cat(x,y,subscripts, lower,upper, length=.2)
}
\arguments{
  \item{x,y}{Data from the call to \code{xyplot}.}
  \item{subscripts}{Variable used to created the juxtaposed panels.}
  \item{lower, upper}{95\% lower and upper bounds of \code{y}.}
  \item{length}{Length of the arrow head lines.}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
library(lattice)
library(effects)
data(Duncan)
Duncan$inc.cat <- cut(Duncan$income, 3)
mod <- lm(prestige~ inc.cat * type + education,
  data=Duncan)
e1 <- effect("inc.cat*type", mod)
update(plot(e1), panel=panel.2cat)
}
