\name{scaleDataFrame}
\alias{scaleDataFrame}
\title{Standardize quantitative variables in a data frame}
\description{This function standardizes quantitative variables in a data frame while leaving the others untouched.  This leaves not only factors, but also binary variables (those with values 0, 1, or NA).}
\usage{
scaleDataFrame(data)
}
\arguments{
  \item{data}{A data frame.}
}
\value{
A data frame with standardized quantitative variables
}
\author{Dave Armstrong}
