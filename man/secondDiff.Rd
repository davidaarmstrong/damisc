% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{secondDiff}
\alias{secondDiff}
\title{Calculate Cross-Derivative and its Variability}
\usage{
secondDiff(
  obj,
  vars,
  data,
  method = c("AME", "MER"),
  vals = NULL,
  typical = NULL
)
}
\arguments{
\item{obj}{An object of class \code{glm} that will be used to find the
cross-derivative.}

\item{vars}{A vector of two variables to be used in calculating the
derivative.}

\item{data}{A data frame.}

\item{method}{Indicate whether you want to use average marginal effects
(AME) or marginal effects at representative values (MER).}

\item{vals}{A named list of length 2 where each element gives the minimum
and maximum values used in the calculation.}

\item{typical}{A named vector of values at which to hold variables constant.}
}
\value{
A list with two elements:

\item{ave}{The average second difference in each iteration of the
bootstrap.} \item{ind}{If \code{type == 'AME'}, \code{ind} is returned with
the second difference and measures of uncertainty for each individual
observation in the original dataset} \item{probs}{If \code{type == 'MER'},
\code{probs} is returned with the full matrix of simulated predicted
probabilities for the four conditions.}
}
\description{
Calculates the cross-derivative required to evaluate interactions in
logistic/probit regression models.
}
\details{
The function calculates the second difference as (Pr(Y=1|x1=max, x2=max) -
Pr(Y=1|x1=min, x2=max)) - (Pr(Y=1|x1=max, x2=min) - Pr(Y=1|x1=min, x2=min)).
The function uses a parametric bootstrap to calculate the sampling
distribution of the second difference.
}
\examples{
data(aclp)
mod <- glm(democ ~ gdpw + popg,
           data = subset(aclp, year == 1980), 
           family=binomial)

secdiff <- secondDiff(mod, c("gdpw", "popg"), data = subset(aclp, year == 1980))
summary(secdiff)

}
\author{
Dave Armstrong
}
