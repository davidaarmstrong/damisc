\name{testNL}
\alias{testNL}
\alias{testNL.glm}
\alias{testNL.lm}
\title{Test Transformations and Polynomials in Non-linear Models}
\description{Tests for model improvements for non-linear transformations and polynomials with Clarke's (2007) distribution-free test for non-nested models. }
\usage{
    \method{testNL}{glm}(obj, var, transPower, polyOrder, plot=FALSE, ...)
    \method{testNL}{lm}(obj, var, transPower, polyOrder, plot=FALSE, ...)
}
\arguments{
    \item{obj}{Object of a supported class in which non-linear functional forms will be tested. }
    \item{var}{String giving name of variable to be tested.}
    \item{transPower}{The power used in the transformation.  For transformations in the range (-0.01, 0.01), the log transformation is used.}
    \item{transPower}{The name value of the power used in the transformation. }
    \item{polyOrder}{The order of the polynomial to be used. }
    \item{plot}{Logical indicating whether the effects should be plotted}
    \item{...}{Currently not implemented.}
}
\details{
Three hypotheses are tested with this function.  The first is whether the original specification is preferred to the power transformation.  The second is whether the original specification is preferred to the polynomial model.  The third is whether the power transformation is preferred to the polynomial model.  All tests are done with the Clarke test. 
}
\value{
A plot or a data frame giving the results of the tests identified above. 
}
\references{
Kevin Clarke.  2007.  "A Simple Distribution-Free Test for
Nonnested Hypotheses."  \emph{Political Analysis} 15(3): 347--363.
}
\author{Dave Armstrong}
