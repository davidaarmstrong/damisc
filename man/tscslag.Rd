\name{tscslag}
\alias{tscslag}
\title{Lag a time-series cross-sectional variables}
\description{Lags (or leads) a variable in a time-series corss-sectional dataset.}
\usage{
	tscslag(dat, x, id, time, lagLength=1)
}
\arguments{
  \item{dat}{A data frame. }
  \item{x}{A string identifying variable to be lagged. }
  \item{id}{A string identifying the name of the cross-sectional identifier. }
  \item{time}{A string identifying the name of the time variable. }
  \item{lagLength}{The length of the lag, use negative values for leading variables.}
}
\value{
A vector giving the lagged values of \code{x}.
}
\author{Dave Armstrong}
