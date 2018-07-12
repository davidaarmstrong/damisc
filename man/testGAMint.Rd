\name{testGAMint}
\alias{testGAMint}
\title{Simulated F-test for Linear Interactions}
\description{Simluates the sampling distribution of the F statistic when comparing a linear intreraction model to a generalized additive model with a smooth over the two variables in the interaction. }
\usage{
testGAMint(m1, m2, data, R=1000, ranCoef=FALSE)
}
\arguments{
  \item{m1}{An object of class \code{gam} estimated with the \code{mgcv} package.  This model should be linear in the interaction of the two x-variables of interest.}
  \item{m2}{An object of class \code{gam} esimtated with the \code{mgcv} package.  This model should contain a smooth interaction.  For two continuous variables, this should be done with \code{te()} unless the variables are measured in the same units (e.g., spatial coordinates) in which case the usual thin-plate regression spline will work.  For categorical moderators, you should use the \code{s(x, by=D0)} and \code{s(x, by=D)} (for a dummy variable moderator, \code{D}, where \code{D0=1} when \code{D=0}.  Remember to include \code{D} as a parametric term in the model as well to account for the intercept difference between the two smooth terms.)}
  \item{data}{Data frame used to estimate both models}
  \item{R}{Number of simulated F values to create.}
  \item{ranCoef}{Logcial indicating whether the coefficients should be treated as fixed or whether they should be drawn from their implied sampling distribution for each iteration of the simulation.}
}
\details{In simple simulations, an F-test of a linear interaction relative to a smooth interaction with a GAM using a nominal .05 type I error rate, has an actual type I error rate of more than double the nominal rate (this tended to be in the low teens).  This function tries to build the F-distribution using simulation.  First, it uses the coefficients from the linear interaction model, multiplies them by the coefficients from the linear interaction model and for each iteration of the simulation, it creates the simulated dependent variable by adding a random error to the linear predictor with the same standard deviation as the residual standard deviation from the linear interaction model.  All of that is to say that this model has all of the same features as the linear interaction model, except that we are certain that this is the right model.  The algorithm then estimates both the linear interaction model and the GAM with a smooth interaction on the original X variables and the new simulated y variable.  The F-test is performed and the F-statistic saved for each iteraction.  The algorithm then calculates the probability of being to the right of the observed F-statistic in the simulated F-distribution.}
\value{
    \item{obsF}{The observed F-statistic from the test on the original models.}
    \item{Fdist}{The \code{R} different F-statistics calculated at each iteration of the simulation.}
}
\author{Dave Armstrong}
