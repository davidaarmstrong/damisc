\name{NKnotsTest}
\alias{NKnotsTest}
\title{Test of functional form assumption using B-splines}
\description{Estimate hypothesis test of lower- and higher-order non-linear relationships against an assumed target relationship.}
\usage{
    NKnotsTest(form, var, data, targetdf=1, degree=3, min.knots=1, 
       max.knots=10, adjust="none")
}
\arguments{
  \item{form}{A formula detailing the model for which smoothing is to be evaluated.}
  \item{var}{A character string identifying the variable for which smoothing is to be evaluated.}
  \item{data}{Data frame providing values of all variables in \code{form}.}
  \item{targetdf}{The assumed degrees of freedom against which the tests will be conducted.}
  \item{degree}{Degree of polynomial in B-spline basis functions.}
  \item{min.knots}{Minimum number of internal B-spline knots to be evaluated.}
  \item{max.knots}{Maximum number of internal B-spline knots to be evaluated.}
  \item{adjust}{Method by which p-values will be adjusted (see \code{\link{p.adjust}})}
}
\value{
    A matrix with the following columns: 
    \item{F}{F statistics of test of candidate models against target model}
    \item{DF1}{Numerator DF from F-test}
    \item{DF2}{Denominator DF from F-test}
    \item{p(F)}{p-value from the F-test}
    \item{Clarke}{Test statistic from the Clarke test}
    \item{Pr(Better)}{The Clarke statistic divided by the number of observations}
    \item{p(Clarke)}{p-value from the Clarke test.  (T) means that the significant p-value is in favor of the Target model and (C) means the significant p-value is in favor of the candidate (alternative) model.}
    \item{Delta_AIC}{AIC(candidate model) - AIC(target model)}
    \item{Delta_AICc}{AICc(candidate model) - AICc(target model)}
    \item{Delta_BIC}{BIC(candidate model) - BIC(target model)}
}
\author{Dave Armstrong}
\examples{
library(car)
NKnotsTest(prestige ~ education + type, var="income", data=na.omit(Prestige), targetdf=3)
}
