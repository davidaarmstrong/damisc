\name{intEff}
\alias{intEff}
\title{Functions for Estimating Interaction Effects in Logit and Probit Models}
\description{Norton and Ai (2003) and Norton, Wang and Ai (2004) discuss methods for calculating the appropriate marginal effects for interactions in binary logit/probit models.  These functions are direct translations of the Norton, Wang and Ai (2004) Stata code.}
\usage{
	intEff(obj, vars, data)
}
\arguments{
  \item{obj}{A binary logit or probit model estimated with \code{glm}.}
  \item{vars}{A vector of the two variables involved in the interaction.}
  \item{data}{A data frame used in the call to \code{obj}.}
}
\value{
A list is returned with two elements - \code{byobs} and \code{atment}.  The \code{byobs} result gives the interaction effect evaluated at each observation.  The \code{atmean} element has the marginal effect evaluated at the mean.  Each eleement contains an element \code{int} which is a data frame with the following variable: 
\item{int_eff}{The correctly calucalted marginal effect.}
\item{linear}{The incorrectly calculated marginal effect following the linear model analogy.}
\item{phat}{Predicted Pr(Y=1|X).}
\item{se_int_eff}{Standard error of \code{int_eff}.}
\item{zstat}{The interaction effect divided by its standard error}

The \code{X} element of each returned result is the X-matrix used to generate the result. 
}
\references{
Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  Computing Interaction Effects and Standard Errors in Logit and Probit Models.  The Stata Journal 4(2): 154-167.\cr

Ai, Chunrong and Edward C. Norton.  2003.  Interaction Terms in Logit and Probit Models.  Economics Letters 80(1): 123-129.

Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  inteff: Computing Interaction Effects and Standard Errors in Logit and Probit Models, Stata Code. 
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
data(france)
mod <- glm(voteleft ~ age*lrself + retnat + male, data=france, family=binomial)
out <- intEff(obj=mod, vars=c("age", "lrself"), data=france)
plot(out$phat, out$int_eff, xlab="Predicted Pr(Y=1|X)", 
	ylab = "Interaction Effect")
ag <- aggregate(out$linear, list(out$phat), mean)
lines(ag[,1], ag[,2], lty=2, col="red", lwd=2)
legend("topright", c("Correct Marginal Effect", "Linear Marginal Effect"), 
	pch=c(1, NA), lty=c(NA, 2), col=c("black", "red"), lwd=c(NA, 2), inset=.01)
}