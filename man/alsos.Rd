% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{alsos}
\alias{alsos}
\title{Alternating Least Squares Optimal Scaling}
\usage{
alsos(
  os_form,
  raw_form = ~1,
  data,
  scale_dv = FALSE,
  maxit = 30,
  level = 2,
  process = 1,
  starts = NULL,
  ...
)
}
\arguments{
\item{os_form}{A two-sided formula including the independent variables to 
be scaled on the left-hand side.  Optionally, the dependent variable can 
also be scaled.}

\item{raw_form}{A right-sided formula with covariates that will not be 
scaled.}

\item{data}{A data frame.}

\item{scale_dv}{Logical indicating whether the dependent variable should 
be optimally scaled.}

\item{maxit}{Maximum number of iterations of the optimal scaling algorithm.}

\item{level}{Measurement level of the dependent variable 1=Nominal,
2=Ordinal}

\item{process}{Nature of the measurement process: 1=discrete, 2=continuous.
Basically identifies whether tied observations will continue to be tied in
the optimally scaled variale (1) or whether the algorithm can untie the
points (2) subject to the overall measurement constraints in the model.}

\item{starts}{Optional starting values for the optimal scaling algorithm.}

\item{...}{Other arguments to be passed down to \code{lm}.}
}
\value{
A list with the following elements:

\item{result}{The result of the optimal scaling process}

\item{data}{The original data frame with additional columns adding the
optimally scaled DV}

\item{iterations}{The iteration history of the algorithm}

\item{form}{Original formula}
}
\description{
Estimates the Alternating Least Squares Optimal Scaling (ALSOS) solution for
qualitative variables.
}
\examples{
data(france)
res <- alsos(lrself ~ 1, ~retnat + male + age, data=france, scale_dv=TRUE)
summary(res$result)

}
\references{
Jacoby, William G.  1999.  \sQuote{Levels of Measurement and Political
Research: An Optimistic View} American Journal of Political Science 43(1):
271-301.

Young, Forrest.  1981.  \sQuote{Quantitative Analysis of Qualitative Data}
Psychometrika, 46: 357-388.

Young, Forrest, Jan de Leeuw and Yoshio Takane.  1976.  \sQuote{Regression
with Qualitative and Quantitative Variables: An Alternating Least Squares
Method with Optimal Scaling Features} Psychometrika, 41:502-529.
}
\author{
Dave Armstrong and Bill Jacoby
}
