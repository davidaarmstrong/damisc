\name{glmChange2}
\alias{glmChange2}
\title{Maximal First Differences for Generalized Linear Models}
\description{For objects of class \code{glm}, it calculates the change in predicted responses, for maximal discrete changes in all covariates holding all other variables constant at typical values.}
\usage{
glmChange2(obj, varname, data, change=c("unit", "sd"), R=1500)
}
\arguments{
  \item{obj}{A model object of class \code{glm}.}
  \item{varname}{Character string giving the variable name for which average effects are to be calculated.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{change}{A string indicating the difference in predictor values to calculate the discrete change.  \code{sd} gives plus and minus one-half standard deviation change around the median and \code{unit} gives a plus and minus one-half unit change around the median.}
  \item{R}{Number of simulations to perform.}
}
\details{
The function calculates the average effect discrete changes in the covariates, for objects of class \code{glm}.  This function works with polynomials specified with the \code{poly} function.  
}
\value{
A vector of values giving the average and 95 percent confidence bounds
}
\author{Dave Armstrong}
\examples{
data(france)
left.mod <- glm(voteleft ~ male + age + retnat + 
	poly(lrself, 2), data=france, family=binomial)
typical.france <- data.frame(
	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
	age = 35
	)
glmChange2(left.mod, "age", data=france, "sd")
}