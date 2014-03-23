\name{aveEffPlot}
\alias{aveEffPlot}
\title{Average Effect Plot for Generalized Linear Models}
\description{For objects of class \code{glm}, it calculates the change the average predicted probability for a hypothetical candidate set of values of a covariate.}
\usage{
aveEffPlot(obj, varname, data, R=1500, nvals=25, plot=TRUE,...) 
}
\arguments{
  \item{obj}{A model object of class \code{glm}.}
  \item{varname}{Character string giving the variable name for which average effects are to be calculated.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{R}{Number of simulations to perform.}
  \item{nvals}{Number of evaluation points at which the average probability will be calculated.}
  \item{plot}{Logical indicating whether plot should be returned, or just data (if \code{FALSE}).}
  \item{...}{Other arguments to be passed down to \code{xyplot}.}
}
\details{
The function plots the average effect of a model covariate, for objects of class \code{glm}.  The function does not work with \code{poly} unless the coefficients are provided as arguments to the command in the model (see example below). 
}
\value{
A plot or a data frame
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
data(france)
p <- poly(france$lrself, 2)
left.mod <- glm(voteleft ~ male + age + retnat + 
	poly(lrself, 2, coefs=attr(p, "coefs")), data=france, family=binomial)
aveEffPlot(left.mod, "age", data=france, plot=FALSE)
}