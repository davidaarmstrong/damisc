\name{probgroup}
\alias{probgroup}
\alias{probgroup.polr}
\alias{probgroup.multinom}
\title{Plog Probabilities by Group}
\description{Plots predicted probabilities by value of the dependent variable for proportional odds logistic regression and multinomial logistic regression models}
\usage{
    \method{probgroup}{polr}(obj, ...)
    \method{probgroup}{multinom}(obj, ...)
}
\arguments{
    \item{obj}{Object of class \code{polr} or \code{multinom} where appropraite. }
    \item{...}{Currently not implemented.}
}
\details{
Plots the predicted probabilities by value of the dependent variable.  Each panel only includes the observations that had the identified value of the dependent variable. 
}
\value{
A plot. 
}
\author{Dave Armstrong}
