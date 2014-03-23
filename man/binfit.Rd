\name{binfit}
\alias{binfit}
\title{Scalar Measures of Fit for Binary Variable Models}
\description{Calculates scalar measures of fit for models with binary dependent variables along the lines described in Long (1997) and Long and Freese (2005).}
\usage{
	binfit(mod)
}
\arguments{
  \item{mod}{A model of class \code{glm} with \code{family=binomial}. }
}
\details{
\code{binfit} calculates scalar measures of fit (many of which are pseudo-R-squared measures) to describe how well a model fits data with a binary dependent variable.  
}
\value{
A named vector of scalar measures of fit
}
\references{
Long, J.S.  1997.  Regression Models for Categorical and Limited Dependent Variables.  Thousand Oaks, CA: Sage.

Long, J.S. and J. Freese.  2005.  Regression Models for Categorical Outcomes Using Stata, 2nd ed.  College Station, TX: Stata Press. 
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
data(france)
left.mod <- glm(voteleft ~ male + age + retnat + 
	poly(lrself, 2), data=france, family=binomial)
binfit(left.mod)
}