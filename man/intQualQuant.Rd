\name{intQualQuant}
\alias{intQualQuant}
\title{Predictions for Factor-Numeric Interactions in Linear Models}
\description{This function works on linear models with a single interaction between a continuous (numeric) variable and a factor.  The output is a data frame that gives the predicted effect of moving from each category to each other category of the factor over the range of values of the continuous conditioning variable.}
\usage{
	intQualQuant(obj, vars, level = .95, labs=NULL, 
		n=10, onlySig=FALSE, plot=FALSE, 
		vals = NULL, rug=TRUE)
}
\arguments{
  \item{obj}{An object of class \code{lm}.}
  \item{vars}{A vector of two variable names giving the two quantitative variables involved in the interaction.  These variables must be involved in one, and only one, interaction. }
  \item{level}{Confidence level desired for lower and upper bounds of confidence interval.}
  \item{labs}{An optional vector of labels that will be used to identify the effects, if \code{NULL}, the factor levels will be used.}
  \item{n}{Number of values of the conditioning variable to use.}
  \item{onlySig}{Logical indicating whether only contrasts with significant differences should be returned.  Significance is determined to exist if the largest lower bound is greater than zero or the smallest upper bound is smaller than zero.}
  \item{plot}{Logical indicating whether a plot should be made (if \code{TRUE}) or the data should be returned (if \code{FALSE}).}	
  \item{vals}{A vector of values at which the continuous variable will be held constant.  If \code{NULL}, a sequence of length \code{n} across the variable's range will be used.}
  \item{rug}{Logical indicating whether rug plots should be plotted in the panels.}
}
\value{
A data frame with the following values: 
	\item{fit}{The expected difference between the two factor levels at the specified value of the conditioning variable.}
	\item{se.fit}{The standard error of the expected differences. }
	\item{x}{The value of the continuous conditioning variable} 
	\item{contrast}{A factor giving the two values of the factor being evaluated.}
	\item{lower}{The lower 95\% confidence interval for \code{fit}}
	\item{upper}{The upper 95\% confidence interval for \code{fit}}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
library(car)
data(Ornstein)
mod <- lm(interlocks ~ log(assets)*nation, data=Ornstein)
mod.out <- intQualQuant(mod, c("log(assets)", "nation"), 
	n=25, plot=FALSE)
}

