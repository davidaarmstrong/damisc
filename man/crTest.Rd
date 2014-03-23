\name{crTest}
\alias{crTest}
\title{Test of linearity for Component + Residual Plots}
\description{This function estimates a linear model and a loess model on the component-plus-residual plot (i.e., a partial residual plot) for each quantitative variable in the model.  The residual sums of squares for each are used to calculate an F-test for each quantitative variable.  }
\usage{
crTest(model, adjust.method="none", ...)
}
\arguments{
  \item{model}{A model object of class \code{lm}}
  \item{adjust.method}{Adjustment method for multiple-testing procedure, using \code{p.adjust} from \code{stats}.}
  \item{...}{Other arguments to be passed down to the call to \code{loess}. }
}
\value{
	A matrix with the following columns for each variable: 
	\item{RSSp}{Residual sum-of-squares for the parametric (linear) model.}
	\item{RSSnp}{Residual sum-of-squares for the non-parametric (loess) model.}
	\item{DFnum}{Numerator degrees of freedom for the F-test: tr(S)-(k+1).}
	\item{DFdenom}{Denominator degrees of freedom for the F-test: n-tr(S)}
	\item{F}{F-statistic}
	\item{p}{p-value, potentially adjusted for multiple comparisons.}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
mod <- lm(prestige ~ income + education + women, data=Prestige)
crTest(mod)
}
