

#' Example data for btscs function
#' 
#' A subset of data from Alvarez et. al. (1996).
#' 
#' 
#' @name aclp
#' @docType data
#' @format A data frame with 4126 observations on the following 7 variables.
#' \describe{ \item{cname}{Country name} \item{country}{Numeric
#' country identifier} \item{year}{Year of observation}
#' \item{reg}{A dichotomous variable coded 1 for dictatorship, 0 for
#' democracy} \item{gdpw}{GDP/worker, 1985 prices}
#' \item{popg}{Population growth} \item{democ}{A dichotomous
#' variable coded 1 for democracy, 0 for dictatorship, (1-\code{reg})} }
#' @references Alvarez, M., J.A. Cheibub, F. Limongi and A. Przeworski. 1996.
#' Classifying political regimes. Studies in Comparative International
#' Development 31 (Summer): 1-37.
#' @keywords datasets
NULL





#' Dave Armstrong's Miscellaneous Functions
#' 
#' Functions to aid in the presentation of linear model results
#' 
#' \tabular{ll}{ Package: \tab DAMisc\cr Type: \tab Package\cr Version: \tab
#' 1.4-8\cr Date: \tab 2018-05-16\cr License: \tab GPL (>=2)\cr LazyLoad: \tab
#' yes\cr } These are functions that help present linear model results.
#' Largely, the represent alternatives in presentation to other R packages.
#' For example, the \code{factorplot} function offers an alternative to David
#' Firth's \code{qvcalc} package.  This function calculates and presents exact
#' variances of all simple contrasts.  Both \code{DAintfun} and
#' \code{DAintfun2} are alternative ways of presenting interactions between two
#' continuous variables.  \code{DAintfun2} gives results in line with the
#' suggestions in Brambor, Clark and Golder (2006).
#' 
#' @name DAMisc-package
#' @aliases DAMisc-package DAMisc
#' @docType package
#' @author Dave Armstrong\cr Maintainer: Dave Armstrong
#' <davearmstrong.ps@@gmail.com>
#' @references Armstrong, D.A.  (2013) factorplot: Improving Presentation of
#' Simple Contrasts in GLMs.  The R Journal.  5, 4-15. Brambor, T., W.R. Clark
#' and M. Golder.  (2006) Understanding Interaction Models: Improving Empirical
#' Analyses.  Political Analysis 14, 63-82.\cr Berrym, W., M. Golder and D.
#' Milton.  (2012) Improving Tests of Theories Positing Interaction.  Journal
#' of Politics 74, 653-671.
#' 
#' @importFrom glue glue
#' @importFrom clarkeTest clarke_test
#' @importFrom grDevices col2rgb dev.off pdf png
#' postscript rgb setEPS
#' @importFrom graphics abline axis lines par persp
#' plot points polygon box segments strwidth
#' @importFrom stats AIC BIC alias anova as.formula
#' coef coefficients contrasts deviance dnorm
#' family fitted formula lm loess logLik
#' median model.frame model.matrix model.response
#' na.omit p.adjust pcauchy pchisq pf plogis
#' pnorm predict pt qt quantile relevel
#' residuals residuals.glm sd terms update
#' var vcov gaussian glm aggregate cooks.distance 
#' cor hatvalues naprint nlminb pbinom printCoefmat
#' qnorm rnorm binomial
#' @importFrom boot cv.glm
#' @importFrom coda as.mcmc
#' @importFrom latticeExtra useOuterStrips
#' @importFrom AICcmodavg AICc
#' @importFrom fANCOVA loess.as
#' @importFrom VGAM yeo.johnson
#' @importFrom optiscale opscale
#' @importFrom utils combn
#' @importFrom VGAM pgumbel
#' @importFrom effects Effect effect
#' @importFrom rstan stan
#' @importFrom MASS mvrnorm polr kde2d glm.nb
#' @importFrom nnet multinom
#' @importFrom splines bs
#' @importFrom pscl zeroinfl
#' @importFrom xtable xtable print.xtable
#' @importFrom ggplot2 ggplot geom_point aes geom_segment 
#' theme_bw labs ggtitle geom_hline geom_ribbon geom_line
#' theme aes_string
#' @importFrom car Anova
#' @importFrom gdata trim
#' @importFrom grid gpar grid.segments unit
#' @importFrom lattice histogram packet.number panel.abline panel.arrows 
#' panel.lines panel.points panel.polygon panel.rug panel.segments 
#' panel.superpose simpleKey trellis.focus trellis.par.get 
#' trellis.unfocus xyplot
NULL





#' Example data for factorplot function
#' 
#' A subset of data from the 1994 Eurobarometer for France
#' 
#' 
#' @name france
#' @docType data
#' @format A data frame with 542 observations on the following 5 variables.
#' \describe{ \item{lrself}{respondent's left-right self-placement on a
#' 1(left)-10(right) scale} \item{male}{a dummy variable coded 1 for
#' males and 0 for females} \item{age}{respondent's age}
#' \item{vote}{a factor indicating vote choice with levels PCF, PS,
#' Green, RPR and UDF} \item{retnat}{a factor indicating the
#' respondent's retrospective national economic evaluation with levels Better,
#' Same and Worse} \item{voteleft}{a dichotomous variable where 1
#' indicates a vote for a left party, 0 otherwise} }
#' @references Reif, Karlheinz and Eric Marlier.  1997. \emph{Euro-barometer
#' 42.0: The First Year of the New European Union, November-December 1994}.
#' Inter-university Consortium for Political and Social Research (ICPSR)
#' [distributor].
#' @keywords datasets
NULL





#' Example Data for DAintfun
#' 
#' Data to execute example code for \code{DAintfun}
#' 
#' These are randomly generated data to highlight the functionality of
#' \code{DAintfun}
#' 
#' @name InteractionEx
#' @docType data
#' @format A data frame with 500 observations on the following 4 variables.
#' \describe{ \item{y}{a numeric vector} \item{x1}{a numeric
#' vector} \item{x2}{a numeric vector} \item{z}{a numeric
#' vector} }
#' @keywords datasets
NULL



