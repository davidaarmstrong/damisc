\name{summary.polr}
\alias{summary.polr}
\alias{print.summary.polr}
\title{Alternative Summary for Ordered Regression Models}
\description{This alternative summary method prints p-values for models estimated with \code{polr}.}
\usage{
\method{summary}{polr}{object, digits = max(3, .Options$digits - 3), correlation = FALSE, ...}
}
\arguments{
  \item{object}{An object of class \code{polr}.}
  \item{digits}{Number of digits to print. }
  \item{clorrelation}{Whether to print the correlation matrix.}
  \item{...}{Other arguments to be passed down}
}
\value{
A model summary
}
\author{Dave Armstrong}
