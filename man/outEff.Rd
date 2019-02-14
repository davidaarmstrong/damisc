\name{outEff}
\alias{outEff}
\title{Plot Effects of Removing Outliers}
\description{Plots the effect of a variable sequentially removing outlying observations.}
\usage{
outEff(obj, var, data, stat =c("cooksD", "hat", "deviance", "pearson"), nOut = 10, whichOut=NULL, cumulative=FALSE)
}
\arguments{
  \item{obj}{An object of class \code{glm} that will be used to plot the outlier removed lines. }
  \item{var}{A character string giving the name of the variable to be used. }
  \item{data}{A data frame.}
  \item{stat}{Which statistic to use to evaluate \sQuote{outlyingness}. }
  \item{nOut}{Number of outliers to be removed.}
  \item{whichOut}{If not \code{NULL}, a vector of observation numbers to be removed manually, rather than using \code{stat}.}
  \item{cumulative}{Logical indicating whether the outliers should be removed cumulatively, or each one in turn, replacing the other outlying observations.}
}
\author{Dave Armstrong}
