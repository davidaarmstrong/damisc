\name{DAintfun}
\alias{DAintfun}
\title{Surface Plots for Two-Way Interactions}
\description{Makes surface plots to display interactions between two continuous variables}
\usage{
DAintfun(obj, varnames, theta = 45, phi = 10, xlab=NULL, ylab=NULL, zlab=NULL,...)
}
\arguments{
  \item{obj}{A model object of class \code{lm}}
  \item{varnames}{A two-element character vector where each element is the name of a variable involved in a two-way interaction.}
  \item{theta}{Angle defining the azimuthal viewing direction to be passed to \code{persp}}
  \item{phi}{Angle defining the colatitude viewing direction to be passed to \code{persp}}
  \item{xlab}{Optional label to put on the x-axis, otherwise if \code{NULL}, it will take the first element of \code{varnames}}
  \item{ylab}{Optional label to put on the y-axis, otherwise if \code{NULL}, it will take the second element of \code{varnames}}
  \item{zlab}{Optional label to put on the z-axis, otherwise if \code{NULL}, it will be \sQuote{Predictions}}
  \item{...}{Other arguments to be passed down to the initial call to \code{persp}}
}
\details{
This function makes a surface plot of an interaction between two continuous covariates.  If the model is \cr
\deqn{y_{i} = b_{0} + b_{1}x_{i1} + b_{2}x_{i2} + b_{3}x_{i1}\times x_{i2} + \ldots + e_{i},}{y = b[0] + b[1]x1 + b[2]x2 + b[3]x1*x2 + ... + e[i],}\cr
this function plots \eqn{b_{1}x_{i1} + b_{2}x_{i2} + b_{3}x_{i1}\times x_{i2}}{b[1]x1 + b[2]x2 + b[3]x1*x2} for values over the range of \eqn{X_{1}}{X1} and \eqn{X_{2}}{X2}.  The highest 75\%, 50\% and 25\% of the bivariate density of \eqn{X_{1}}{X1} and \eqn{X_{2}}{X2} (as calculated by \code{sm.density} from the \code{sm} package) are colored in with colors of increasing gray-scale. 
}
\value{
\item{x1}{Values of the first element of \code{varnames} used to make predictions.}
\item{x2}{Values of the second element of \code{varnames} used to make predictions.}
\item{pred}{The predictions based on the values \code{x1} and \code{x2}.}
\item{graph}{A graph is produced, but no other information is returned.}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
data(InteractionEx)
mod <- lm(y ~ x1*x2 + z, data=InteractionEx)
DAintfun(mod, c("x1", "x2"))

## Make interactive with:
# mypanel <- function(panel){
# 	DAintfun(mod, c("x1", "x2"), theta=panel$theta, phi=panel$phi)
# 	panel}
# panel <- rp.control(theta=0, phi=25)
# rp.slider(panel, theta, -360, 360, mypanel, showvalue=TRUE)
# rp.slider(panel, phi, 0, 90, mypanel, showvalue=TRUE)
}