\name{BGMtest}
\alias{BGMtest}
\title{Tests the five Berry, Golder and Milton (2012) Interactive Hypothesis}
\description{This function tests the five hypotheses that Berry, Golder and Milton identify as important when two quantitative variables are interacted in a linear model.}
\usage{
BGMtest(obj, vars, digits = 3, level = 0.05, two.sided=T)
}
\arguments{
  \item{obj}{An object of class \code{lm}.}
  \item{vars}{A vector of two variable names giving the two quantitative variables involved in the interaction.  These variables must be involved in one, and only one, interaction. }
  \item{digits}{Number of digits to be printed in the summary.}
  \item{level}{Type I error rate for the tests.}
  \item{two.sided}{Logical indicating whether the tests should be two-sided (if \code{TRUE}, the default) or one-sided (if \code{FALSE}).}
}
\value{
A matrix giving five t-tests. 
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
data(Duncan)
mod <- lm(prestige ~ income*education + type, data=Duncan)
BGMtest(mod, c("income", "education"))
}