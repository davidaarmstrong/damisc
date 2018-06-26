\name{combTest}
\alias{combTest}
\title{Test for Combining Categories in Multinomial Logistic Regression Models.}
\description{Tests the null hypothesis that categories can be combined in Multinomial Logistic Regression Models}
\usage{
   combTest(obj)
}
\arguments{
  \item{obj}{An object of class \code{multinom}.}
}
\value{
A matrix of test statistics and p-values. 
}
\author{Dave Armstrong}
\examples{
library(nnet)
data(france)
mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
combTest(mnl.mod)

}
