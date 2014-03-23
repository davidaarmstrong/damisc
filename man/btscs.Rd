\name{btscs}
\alias{btscs}
\title{Generate Spells for Binary Variables}
\description{Beck et. al. (1998) identified that binary time-series cross-section data are discrete-time duration data and time dependence can be modeled in a logistic regression by including a flexible function (e.g., cubic spline) of time since the last event as a covariate.  This function creates the variable identifying time since last event.}
\usage{
btscs(data, event, tvar, csunit, pad.ts=FALSE)
}
\arguments{
  \item{data}{A data frame.}
  \item{event}{Character string giving the name of the dichotomous variable identifying the event (where an event is coded 1 and the absence of an event is coded 0).}
  \item{tvar}{Character string giving the name of the time variable.}
  \item{csunit}{Character string giving the name of the cross-sectional unit.}
  \item{pad.ts}{Logical indicating whether the time-series should be filled in, when panels are unbalanced.}
}
\value{
The original data frame with one additional variable.  The \code{spell} variable identifies the number of observed periods since the last event. 
}
\references{
Alvarez, M., J.A. Cheibub, F. Limongi and A. Przeworski. 1996. Classifying political regimes. Studies in Comparative International Development 31 (Summer): 1-37.

Beck, N.. J. Katz and R. Tucker. 1998. Beyond Ordinary Logit: Taking Time Seriously in Binary-Time-Series-Cross-Section Models. American Journal of Political Science 42(4): 1260-1288.
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(splines)
## Data from Alvarez et. al. (1996)
data(aclp)
newdat <- btscs(aclp, "democ", "year", "country")	

# Estimate Model with and without spell
full.mod <- glm(democ ~ log(gdpw) + popg + bs(spell, df=4), data=newdat, family=binomial)	
restricted.mod <- glm(democ ~ log(gdpw) + popg, data=newdat, family=binomial)	
	
# Incremental F-test of time dependence
anova(restricted.mod, full.mod, test='Chisq')	
}