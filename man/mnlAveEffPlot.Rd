\name{mnlAveEffPlot}
\alias{mnlAveEffPlot}
\title{Average Effects Plot for Multinomial Logistic Regression}
\description{Produces a plot of average effects for one variable while holding the others constant at observed values.}
\usage{
   mnlAveEffPlot(obj, varname, data, R = 1500, nvals = 25, plot = TRUE, ...)
}
\arguments{
  \item{obj}{An object of class \code{multinom}.}
  \item{varname}{A string indicating the variable for which the plot is desired.}
  \item{data}{The data used to estimate \code{obj}.}
  \item{R}{Number of simulations used to generate confidence bounds.}
  \item{nvals}{Number of evaluation points for the predicted probabilities.}
  \item{plot}{Logical indicating whether a plot should be produced (if \code{TRUE}) or numerical results should be returned (if \code{FALSE}).}
  \item{\dots}{Other arguments to be passed down to \code{xyplot}.}
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
\author{Dave Armstrong}

\examples{
library(nnet)
data(france)
mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
\dontrun{mnlAveEffPlot(mnl.mod, "lrself", data=france)}
}
