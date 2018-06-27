\name{yj_trans}
\alias{yj_trans}
\title{Optimizing Yeo-Johnson Transformation}
\description{Uses \code{nlminb} to find the optimal Yeo-Johnson transformation parameters conditional on a parametric model specification.}
\usage{
	yj_trans(form, data, trans.vars, round.digits=3, ...)
}
\arguments{
  \item{form}{A formula with a dependent variable that will be optimally scaled}

  \item{data}{A data frame.}

  \item{trans.vars}{A character string identifying the variables that should be transformed}

  \item{round.digits}{Number of digits to round the transformation parameters.}

  \item{...}{Other arguments to be passed down to \code{lm}.}
}
\value{
A linear model object that was estimated on the optimally transformed variables. 
}
\author{Dave Armstrong}
