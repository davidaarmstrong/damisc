\name{simPredpolr}
\alias{simPredpolr}
\title{Calculate Predictions for Proportional Odds Logistic Regression}
\description{Calculates predicted probabilities from models of class \code{polr} from a model object and a vector of coefficient values.  This is an auxiliary function used in \code{pre} if \code{sim=TRUE}.}
\usage{
simPredpolr(object, coefs, n.coef) 
}
\arguments{
  \item{object}{An object of class \code{polr}.}
  \item{n.coef}{Number of coefficients (minus intercepts) for the \code{polr} model.}
  \item{coefs}{A vector of coefficients where elements 1 to \code{n.coef} give model coefficients and elements \code{n.coef}+1 to k have intercepts.}
}
\value{
An n x m-category matrix of predicted probabilities
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
