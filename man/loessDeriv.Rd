\name{loessDeriv}
\alias{loessDeriv}
\title{Estimate Derivatives of LOESS Curve.}
\description{Estimates the first derivatives of the LOESS curve.}
\usage{
	loessDeriv(obj, delta=.00001)
}
\arguments{
  \item{obj}{An object of class \code{loess}.}

  \item{delta}{Small change to be induced to estimate derivative. }

}
\value{
A vector of first derivative values evaluated at each original x-value. 
}
\author{Dave Armstrong}
