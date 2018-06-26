\name{ordChange2}
\alias{ordChange2}
\title{Average Effects for Proportional Odds Logistic Regression Models}
\description{For objects of class \code{polr}, it calculates the average change in predicted probabilities, for discrete changes in a covariate holding all other variables at their observed values.}
\usage{
ordChange2(obj, varnames, data, diffchange = c("sd", "unit"), 
   R=1500)
}
\arguments{
  \item{obj}{A model object of class \code{polr}.}
  \item{varnames}{A vector of strings identifying the variable to be manipulated.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{diffchange}{A string indicating the difference in predictor values to calculate the discrete change. \code{sd} gives plus and minus one-half standard deviation change around the median and \code{unit} gives a plus and minus one-half unit change around the median.}
  \item{R}{Number of simulations.}
}
\details{
The function calculates the changes in predicted probabilities for maximal discrete changes in the covariates for objects of class \code{polr}.  This function works with polynomials specified with the \code{poly} function.  It also works with multiplicative interactions of the covariates by virtue of the fact that it holds all other variables at typical values.  By default, typical values are the median for quantitative variables and the mode for factors. The way the function works with factors is a bit different.  The function identifies the two most different levels of the factor and calculates the change in predictions for a change from the level with the smallest prediction to the level with the largest prediction.  
}
\value{
A list with the following elements:
	\item{diffs}{A matrix of calculated first differences}
	\item{minmax}{A matrix of values that were used to calculate the predicted changes}
	\item{minPred}{A matrix of predicted probabilities when each variable is held at its minimum value, in turn.}
	\item{maxPred}{A matrix of predicted probabilities when each variable is held at its maximum value, in turn.}
}
\author{Dave Armstrong}
\examples{
library(MASS)
data(france)
polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
typical.france <- data.frame(
	age = 35, 
	retnat = factor(1, levels=1:3, labels=levels(france$retnat))
	)
ordChange2(polr.mod, "age", data=france, diffchange="sd")	
}