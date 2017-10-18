\name{searchVarLabels}
\alias{searchVarLabels}
\alias{searchVarLabels.data.frame}
\alias{searchVarLabels.tbl_df}
\title{Search Variable Labels Attribute}
\description{Data imported from SPSS or Stata comes with the variable labels set (if they were set in the original dataset) as one of the dataframe's attributes.  This allows you to search the variable labels and returns the variable column number, name and label for all variables that have partially match the search term either in their labels or names.}
\usage{
searchVarLabels(dat, str)
}
\arguments{
  \item{dat}{a data frame whose variable labels you want to search.}
  \item{str}{string used to search variable labels.}
}
\details{
For an imported Stata dataset, variable labels are in the \code{var.labels} attribute of the dataset and in an SPSS dataset, they are in the \code{variable.labels} attribute.  These are searched, ignoring case, for the desired string
}
\value{
\item{matrix}{A matrix of dimensions n-matches x 2 is returned, where the first column is the column number of the matching variable and the second column is the variable label.  The row names of the matrix are the variable names.}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
