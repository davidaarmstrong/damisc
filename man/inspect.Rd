\name{inspect}
\alias{inspect}
\alias{inspect.data.frame}
\alias{inspect.tbl_df}
\title{Inspect a Variable in a Data Frame}
\description{Shows the variable label, factor levels (i.e., value labels) and frequency distribution for the specified variable.}
\usage{
   inspect(data, x, ...)
}
\arguments{
  \item{data}{A data frame of class \code{data.frame} or \code{tbl_df}.}

  \item{x}{A string identifying the name of the variable to be inspected.}

  \item{...}{Other arguments to be passed down, currently unimplemented.}
}
\value{
A list with a variable label (if present), factor levels/value labels (if present) and a frequency distribution
}
\author{Dave Armstrong}
\examples{
data(france)
inspect(france, "vote")
}
