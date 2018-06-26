\name{intQualQuant}
\alias{intQualQuant}
\title{Predictions for Factor-Numeric Interactions in Linear Models}
\description{This function works on linear models with a single interaction between a continuous (numeric) variable and a factor.  The output is a data frame that gives the predicted effect of moving from each category to each other category of the factor over the range of values of the continuous conditioning variable.}
\usage{
	intQualQuant(obj, vars, level = .95, varcov=NULL, labs=NULL, 
		n=10, onlySig=FALSE, type=c("facs", "slopes"),
		plot=TRUE, vals = NULL, rug=TRUE, ci=TRUE, digits=3,...)
}
\arguments{
  \item{obj}{An object of class \code{lm}.}
  \item{vars}{A vector of two variable names giving the two quantitative variables involved in the interaction.  These variables must be involved in one, and only one, interaction. }
  \item{level}{Confidence level desired for lower and upper bounds of confidence interval.}
  \item{varcov}{A potentially clustered or robust variance-covariance matrix of parameters used to calculate standard errors.  If \code{NULL}, the \code{vcov} function will be used. }
  \item{labs}{An optional vector of labels that will be used to identify the effects, if \code{NULL}, the factor levels will be used.}
  \item{n}{Number of values of the conditioning variable to use.}
  \item{onlySig}{Logical indicating whether only contrasts with significant differences should be returned.  Significance is determined to exist if the largest lower bound is greater than zero or the smallest upper bound is smaller than zero.}
  \item{type}{String indicating whether the conditional partial effect of the factors is plotted (if \sQuote{facs}), or the conditional partial effect of the quantitative variable (if \sQuote{slopes}) is produced.}	
  \item{plot}{Logical indicating whether graphical results (if \code{TRUE}) or numerical results (if \code{FALSE}) are produced.}
  \item{vals}{A vector of values at which the continuous variable will be held constant.  If \code{NULL}, a sequence of length \code{n} across the variable's range will be used.}
  \item{rug}{Logical indicating whether rug plots should be plotted in the panels.}
  \item{ci}{Logical indicating whether confidence bounds should be drawn.}
  \item{digits}{Number indicating how many decimal places to round the numeric output.}
  \item{...}{Other arguments to be passed down to \code{effect} if \code{plot.type} = \sQuote{slopes}.}
}
\value{
For \code{type} = \sQuote{facs} and \code{plot} = \code{FALSE}, a data frame with the following values: 
	\item{fit}{The expected difference between the two factor levels at the specified value of the conditioning variable.}
	\item{se.fit}{The standard error of the expected differences. }
	\item{x}{The value of the continuous conditioning variable} 
	\item{contrast}{A factor giving the two values of the factor being evaluated.}
	\item{lower}{The lower 95\% confidence interval for \code{fit}}
	\item{upper}{The upper 95\% confidence interval for \code{fit}}
For \code{type} = \sQuote{facs} and \code{plot} = \code{TRUE}, a lattice display is returned
For \code{type} = \sQuote{slopes} and \code{plot} = \code{FALSE}, A character matrix with the following columns: 
	\item{B}{The conditional effect of the quantitative variable for each level of the factor.}
	\item{SE(B)}{The standard error of the conditional effect.}
	\item{t-stat}{The t-statistic of the conditional effect.}
	\item{Pr(>|t|)}{The two-sided p-value.}
For \code{type} = \sQuote{slopes} and \code{plot} = \code{TRUE}, a lattice display is returned
}
\author{Dave Armstrong}
\examples{
library(car)
data(Prestige)
Prestige$income <- Prestige$income/1000
mod <- lm(prestige ~ income * type + education, data=Prestige)
intQualQuant(mod, c("income", "type"), n=10, 
	plot.type="none")
intQualQuant(mod, c("income", "type"), n=10, 
	plot.type="facs")
intQualQuant(mod, c("income", "type"), n=10, 
	plot.type="slopes")
}

