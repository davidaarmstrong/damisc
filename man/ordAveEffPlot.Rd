\name{ordAveEffPlot}
\alias{ordAveEffPlot}
\title{Plot Average Effects of Variables in Proportional Odds Logistic Regression}
\description{
For objects of class \code{polr} the function plots the average effect of a single variable holding all other variables at their observed values. 
}
\usage{
ordAveEffPlot(obj, varname, data, R = 1500, nvals = 25, plot = TRUE, ...)
}
\arguments{
  \item{obj}{An object of class \code{polr}}
  \item{varname}{A string providing the name of the variable for which you want the plot to be drawn.}
  \item{data}{Data used to estimate \code{obj}.}
  \item{R}{Number of simulations to generate confidence intervals.}
  \item{nvals}{Number of evaluation points of the function}
  \item{plot}{Logical indicating whether or not the result should be plotted (if \code{TRUE}) or returned to the console (if \code{FALSE}).}
  \item{\dots}{Arguments passed down to the call to \code{xyplot}}
}
\details{
Following the advice of Hanmer and Kalkan (2013) the function calculates the average effect of a variable holding all other variables at observed values and then plots the result. 
}
\value{
Either a plot or a data frame with variables 
   \item{mean}{The average effect (i.e., predicted probability)}
   \item{lower}{The lower 95\% confidence bound}
   \item{upper}{The upper 95\% confidence bound}
   \item{y}{The values of the dependent variable being predicted}
   \item{x}{The values of the independent variable being manipulated}
}
\references{
Hanmer, M.J. and K.O. Kalkan.  2013.  \sQuote{Behind the Curve: Clarifying the Best Approach to Calculating Predicted Probabilities and Marginal Effects from Limited Dependent Variable Models}.  American Journal of Political Science.  57(1): 263-277.
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}

\examples{
library(MASS)
data(france)
polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
\dontrun{ordAveEffPlot(polr.mod, "lrself", data=france)	}
}
