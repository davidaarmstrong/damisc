\name{glmChange}
\alias{glmChange}
\title{Maximal First Differences for Generalized Linear Models}
\description{For objects of class \code{glm}, it calculates the change in predicted responses, for maximal discrete changes in all covariates holding all other variables constant at typical values.}
\usage{
glmChange(obj, data, typical.dat=NULL, diffchange=c("range", "sd", "unit"), 
    sim=FALSE, R=1000)
}
\arguments{
  \item{obj}{A model object of class \code{glm}.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{typical.dat}{Data frame with a single row containing values at which to hold variables constant when calculating first differences.  These values will be passed to \code{predict}, so factors must take on a single value, but have all possible levels as their levels attribute.}
  \item{diffchange}{A string indicating the difference in predictor values to calculate the discrete change.  \code{range} gives the difference between the minimum and maximum, \code{sd} gives plus and minus one-half standard deviation change around the median and \code{unit} gives a plus and minus one-half unit change around the median.}
  \item{sim}{Logical indicating whether simulated confidence bounds on the difference should be calculated and presented.}
  \item{R}{Number of simulations to perform if \code{sim} is \code{TRUE}}
}
\details{
The function calculates the changes in predicted responses for maximal discrete changes in the covariates, for objects of class \code{glm}.  This function works with polynomials specified with the \code{poly} function.  It also works with multiplicative interactions of the covariates by virtue of the fact that it holds all other variables at typical values.  By default, typical values are the median for quantitative variables and the mode for factors. The way the function works with factors is a bit different.  The function identifies the two most different levels of the factor and calculates the change in predictions for a change from the level with the smallest prediction to the level with the largest prediction.  
}
\value{
A list with the following elements:
	\item{diffs}{A matrix of calculated first differences}
	\item{minmax}{A matrix of values that were used to calculate the predicted changes}
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
glmChange(left.mod, data=france, typical.dat=typical.france)
}