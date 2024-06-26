% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{aveEffPlot}
\alias{aveEffPlot}
\title{Average Effect Plot for Generalized Linear Models}
\usage{
aveEffPlot(
  obj,
  varname,
  data,
  R = 1500,
  nvals = 25,
  level = 0.95,
  ciType = c("percent", "normal"),
  return = c("ci", "plot", "sim"),
  ...
)
}
\arguments{
\item{obj}{A model object of class \code{glm}.}

\item{varname}{Character string giving the variable name for which average
effects are to be calculated.}

\item{data}{Data frame used to fit \code{object}.}

\item{R}{Number of simulations to perform.}

\item{nvals}{Number of evaluation points at which the average probability
will be calculated.}

\item{level}{Scalar giving the confidence level of the point-wise confidence 
intervals.}

\item{ciType}{Type of confidence interval to be created.  If \code{"perc"}, a 
percentile interval will be created from the distribution of effects.  If 
\code{"normal"} a normal-theory interval will be calculated using the standard
deviation of the fitted response from the simulation.}

\item{return}{Character string indicating what should be returned.  Multiple 
entries are supported.}

\item{...}{Other arguments to be passed down to \code{xyplot}.}
}
\value{
A plot or a list depending on the value of `plot`.  If: 
\item{TRUE}{The plot gives the
average predicted probability and confidence bounds moving across the values  
of the independent variable of interest}
\item{FALSE}{A list with three elements - `ci` which has the values of the independent
variable (labelled `s`), the predicted probability in `mean` and the lower and upper
confidence bounds in `lower` and `upper` respectively; an empty element named `plot` and
an element named `sim` that contains the actual simulated values that can be recovered if
desired.}
}
\description{
For objects of class \code{glm}, it calculates the change the average
predicted probability (like the one calculated by \code{glmChange2}) for a
hypothetical candidate set of values of a covariate.
}
\details{
The function plots the average effect of a model covariate, for objects of
class \code{glm}.  The function does not work with \code{poly} unless the
coefficients are provided as arguments to the command in the model (see
example below).
}
\examples{

data(france)
p <- poly(france$lrself, 2)
left.mod <- glm(voteleft ~ male + age + retnat + 
	poly(lrself, 2, coefs=attr(p, "coefs")), data=france, family=binomial)
aveEffPlot(left.mod, "age", data=france, plot=FALSE)

}
\author{
Dave Armstrong
}
