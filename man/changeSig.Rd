\name{changeSig}
\alias{changeSig}
\title{Regions of Statistical Significance in Interactions}
\description{Calculates the regions of statistical significance in interaction and identifies the points at which the statistical significance of conditional coefficients changes. }
\usage{
	changeSig(obj, vars, alpha=.05)
}
\arguments{
  \item{obj}{A model of class \code{glm} or class \code{lm}. }

  \item{vars}{A character vector of the names of the two variables involved in the interaction.}

  \item{alpha}{Critical p-value of the test.}
}
\value{
Printed output that identifies the change-points in statistical significance. 
}
\author{Dave Armstrong}
