\name{ziChange}
\alias{ziChange}
\title{Maximal First Differences for Zero-Inflated Models}
\description{Calculates the change in predicted counts or optionally the predicted probability of being in the zero-count group, for maximal discrete changes in all covariates holding all other variables constant at typical values.}
\usage{
ziChange(obj, data, typical.dat=NULL, type="count")
}
\arguments{
  \item{obj}{A model object of class \code{zeroinfl}.}
  \item{data}{Data frame used to fit \code{object}.}
  \item{typical.dat}{Data frame with a single row containing values at which to hold variables constant when calculating first differences.  These values will be passed to \code{predict}, so factors must take on a single value, but have all possible levels as their levels attribute.}
  \item{type}{Character string of either \sQuote{count} (to obtain changes in predicted counts) or \sQuote{zero} (to obtain changes in the predicted probability of membership in the zero group).}
}
\details{
The function calculates the changes in predicted counts, or optionally the predicted probability of being in the zero group, for maximal discrete changes in the covariates.  This function works with polynomials specified with the \code{poly} function.  It also works with multiplicative interactions of the covariates by virtue of the fact that it holds all other variables at typical values.  By default, typical values are the median for quantitative variables and the mode for factors. The way the function works with factors is a bit different.  The function identifies the two most different levels of the factor and calculates the change in predictions for a change from the level with the smallest prediction to the level with the largest prediction.  
}
\value{
A list with the following elements:
	\item{diffs}{A matrix of calculated first differences}
	\item{minmax}{A matrix of values that were used to calculate the predicted changes}
}
\author{Dave Armstrong}
\examples{
library(pscl)
## Example from the help file for zeroinfl in the pscl package
data("bioChemists", package = "pscl")
fm_zinb <- zeroinfl(art ~ fem + mar + kid5 + phd + ment |
    fem + mar + kid5 + phd + ment, data = bioChemists, dist = "negbin")
typical.bioChem <- data.frame(
	kid5 = 2,
	mar = factor(1, levels=1:2, labels=levels(bioChemists$mar))
	)
ziChange(fm_zinb, data=bioChemists, typical.dat=typical.bioChem, type="zero")	
ziChange(fm_zinb, data=bioChemists, typical.dat=typical.bioChem, type="count")	
}