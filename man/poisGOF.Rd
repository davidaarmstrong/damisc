\name{poisGOF}
\alias{poisGOF}
\title{Deviance and Chi-squared Goodness-of-Fit Test for Poisson Models}
\description{Deviance and Chi-squared goodness-of-fit test of the null hypothesis that poisson variance is appropriate to model the conditional dispersion of the data, given a particular model.}
\usage{
poisGOF(obj)
}
\arguments{
  \item{obj}{A model object of class \code{glm} (with \code{family=poisson}).}
}
\value{
A 2x2 data frame with rows representing the different types of statistics (Deviance and Chi-squared) and columns representing the test statistic and p-value. 
}
\references{
Dobson, A. J. (1990) An Introduction to Generalized Linear Models. London: Chapman and Hall.
}
\author{Dave Armstrong}
\examples{
## Example taken from MASS help file for glm, identified to be
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family=poisson())
poisGOF(glm.D93)
}
