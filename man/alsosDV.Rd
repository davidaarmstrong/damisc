\name{alsosDV}
\alias{alsosDV}
\title{Alternating Least Squares Optimal Scaling}
\description{Estimates the Alternating Least Squares Optimal Scaling (ALSOS) solution for qualitative dependent variables.}
\usage{
	alsosDV(form, data, maxit=30, level=2, process=1, starts=NULL,...)
}
\arguments{
  \item{formula}{A formula with a dependent variable that will be optimally scaled}

  \item{data}{A data frame.}

  \item{maxit}{Maximum number of iterations of the optimal scaling algorithm. }

  \item{level}{Measurement level of the dependent variable 1=Nominal, 2-Ordinal}

  \item{process}{Nature of the measurement process: 1=discrete, 2=continuous.  Basically identifies whether tied observations will continue to be tied in the optimally scaled variale (1) or whether the algorithm can untie the points (2) subject to the overall measurement constraints in the model.}
  
  \item{...}{Other arguments to be passed down to \code{lm}.}
}
\details{
\code{alsosDV} estimates the Alternating Least Squares Optimal Scaling solution on the dependent variable.  
}
\value{
A list with the following elements: 

\item{result}{The result of the optimal scaling process}

\item{data}{The original data frame with additional columns adding the optimally scaled DV}

\item{iterations}{The iteration history of the algorithm}

\item{form}{Original formula}
}

\references{

Jacoby, William G. 1999. \sQuote{Levels of Measurement and Political Research: An Optimistic View.} American Journal of Political Science 43(1):271–301.

Young, Forrest.  1981.  \sQuote{Quantitative Analysis of Qualitative Data}, Psychometrika, 46:357-388. 

Young, Forrest, Jan de Leeuw and Yoshio Takane.  1976.  \sQuote{Regression with Qualitative and Quantitative Variables: An Alternating Least Squares Method with Optimal Scaling Features}.  Psychometrika, 41:502-529.

}
\author{Dave Armstrong}
