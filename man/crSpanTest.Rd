\name{crSpanTest}
\alias{crSpanTest}
\title{Test of Span Parameter in linearity for Component + Residual Plots}
\description{This function performs \code{crTest} for a user-defined range of span parameters, optionally allowing for multiple testing corrections in the p-values. }
\usage{
crSpanTest(model, spfromto, n=10, adjust.method="none", 
    adjust.type=c("none", "across", "within", "both"))
}
\arguments{
  \item{model}{A model object of class \code{lm}}
  \item{spfromto}{A vector of two values across which a range of \code{n} span values will be generated and tested.}
  \item{n}{Number of span parameters to test.}
  \item{adjust.method}{Adjustment method for multiple-testing procedure, using \code{p.adjust} from \code{stats}.}
  \item{adjust.type}{String giving the values over which the  multiple testing correction will be performed.  Here, \sQuote{both} refers to a multiple testing correction done over all span parameters and all variables in the model.  \sQuote{within} means the multiple testing correction should be done within each model, but not across the span parameters and \sQuote{across} means that the multiple testing correction should be for each variable across the various span parameters, but not across variables within the same model.  \sQuote{none} refers to a pass-through option of no multiple testing procedure.}
  \item{...}{Other arguments to be passed down to the call to \code{loess}. }
}
\value{
	A list with two elements: 
	\item{x}{Sequence of span values used in testing}
	\item{y}{p-values for each variable for each span parameter}
}
\author{Dave Armstrong}
\examples{
library(car)
mod <- lm(prestige ~ income + education + women, data=Prestige)
tmp <- crSpanTest(mod, c(.1, .9), adjust.method="holm", 
	adjust.type="both")
matplot(tmp$x, tmp$y, type="l")
}
