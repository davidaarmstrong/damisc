\name{makeHypSurv}
\alias{makeHypSurv}
\title{Make Hypothetical Predictions for Survey Data}
\description{Calculates survival probabilities for hypothetical data.}
\usage{
	makeHypSurv(l, obj, ...)
}
\arguments{
  \item{l}{A named list where variable names are the names and values of the variables are the values.  Combinations will be made with \code{expand.grid}.}
  \item{obj}{A model object estimated with \code{survreg}.}
  \item{...}{currently not implemented.}
}
\value{
A data frame.
}
\author{Dave Armstrong}
