\name{effect.logistf}
\alias{effect.logistf}
\title{Plot Effects from Firth Logit}
\description{Plots the effect of a variable in a model estimated with Firth Logit.}
\usage{
effect.logistf(var, obj, data, ...)
}
\arguments{
  \item{var}{A character string giving the name of the variable whose effect is to be generated. }
  \item{obj}{An object of class \code{logistf}. }
  \item{data}{A data frame.}
  \item{...}{Other arguments to be passed down to the \code{\link{Effect}} function.}
}
\details{
The \code{effect.logistf} function calculates the effect (predicted probabilities) of a variable in a Firth logit model estimated with the \code{\link{logistf}} function.  The function estimates the analogous glm.  It then replaces the coefficient vector in that model object with the Frith logit coefficients.  It also puts the variance-covariance matrix from the Firth logit in the model object and uses a custom extractor function in the \code{Effect} function to extract that variance-covariance matrix rather than the one usually extracted with \code{vcov}.  Note that variability and confidence intervals for the effects will not be calculated using profile likelihood as they are in the Firth logit, but will be calculated using the appropriate variance-covariance matrix.  
}

\value{
An object of class \code{eff} that can be used with other functions from the \code{effects} package. 
}
\author{Dave Armstrong}
