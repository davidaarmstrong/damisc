\name{mnlChange2}
\alias{mnlChange2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Average Effects for Multinomial Logistic Regression Models}
\description{
Calculates average effects of a variable in multinomial logistic regression holding all other variables at observed values. 
}
\usage{
mnlChange2(obj, varname, data, change = c("unit", "sd"), R = 1500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class \code{multinom}}
  \item{varname}{A string identifying the variable to be manipulated.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{change}{A string indicating the difference in predictor values to calculate the discrete change. \code{sd} gives plus and minus one-half standard deviation change around the median and \code{unit} gives a plus and minus one-half unit change around the median.}
  \item{R}{Number of simulations.}
}
\value{
A list with elements: 
   \item{mean}{Average effect of the variable for each category of the dependent variable.}
   \item{lower}{Lower 95 percent confidence bound}
   \item{upper}{Upper 95 percent confidence bound}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(nnet)
data(france)
mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
mnlChange2(mnl.mod, "lrself", data=france, )	


}
