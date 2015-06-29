\name{logit_cc}
\alias{logit_cc}
\alias{logit_cd}
\alias{logit_dd}
\alias{probit_cc}
\alias{probit_cd}
\alias{probit_dd}
\title{Functions for Estimating Interaction Effects in Logit and Probit Models}
\description{Norton and Ai (2003) and Norton, Wang and Ai (2004) discuss methods for calculating the appropriate marginal effects for interactions in binary logit/probit models.  These functions are direct translations of the Norton, Wang and Ai (2004) Stata code.  These functions are not intended to be called by the user directly, rather they are called as needed by \code{intEff}.}
\usage{
   logit_cc(obj=obj, int.var=int.var, vars=vars, b=b, X=X) 
   logit_cd(obj=obj, int.var=int.var, vars=vars, b=b, X=X) 
   logit_dd(obj=obj, int.var=int.var, vars=vars, b=b, X=X)
   probit_cc(obj=obj, int.var=int.var, vars=vars, b=b, X=X) 
   probit_cd(obj=obj, int.var=int.var, vars=vars, b=b, X=X) 
   probit_dd(obj=obj, int.var=int.var, vars=vars, b=b, X=X)
}
\arguments{
  \item{obj}{A binary logit or probit model estimated with \code{glm}.}
  \item{int.var}{The name of the interaction variable.}
  \item{vars}{A vector of the two variables involved in the interaction.}
  \item{b}{Coefficients from the \code{glm} object.}
  \item{X}{Model matrix from the \code{glm} object.}
}
\value{
A data frame with the following variable: 
\item{int_eff}{The correctly calucalted marginal effect.}
\item{linear}{The incorrectly calculated marginal effect following the linear model analogy.}
\item{phat}{Predicted Pr(Y=1|X).}
\item{se_int_eff}{Standard error of \code{int_eff}.}
\item{zstat}{The interaction effect divided by its standard error}
}
\references{
Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  Computing Interaction Effects and Standard Errors in Logit and Probit Models.  The Stata Journal 4(2): 154-167.\cr

Ai, Chunrong and Edward C. Norton.  2003.  Interaction Terms in Logit and Probit Models.  Economics Letters 80(1): 123-129.

Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  inteff: Computing Interaction Effects and Standard Errors in Logit and Probit Models, Stata Code.
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
