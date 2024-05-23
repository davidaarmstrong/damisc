#' Bayesian Alternating Least Squares Optimal Scaling
#' 
#' Estimates a Bayesian analog to the the Alternating Least Squares Optimal
#' Scaling (ALSOS) solution for qualitative dependent variables.
#' 
#' \code{balsos} estimates a Bayesian analog to the Alternating Least Squares
#' Optimal Scaling solution on the dependent variable.  This permits testing
#' linearity assumptions on the original scale of the dependent variable.
#' 
#' @param formula A formula with a dependent variable that will be optimally
#' scaled
#' @param data A data frame.
#' @param iter Number of samples for the MCMC sampler.
#' @param chains Number of parallel chains to be run.
#' @param alg Algorithm used to do sampling.  See \code{stan} for more
#' details.
#' @param ... Other arguments to be passed down to \code{stanfit}.
#' @returns A list with the following elements:
#' 
#' \item{fit}{The fitted stan output}
#' 
#' \item{y}{The dependent variable values used in the regression. }
#' 
#' \item{X}{The design matrix for the regression}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @references
#' 
#' Jacoby, William G.  1999.  \sQuote{Levels of Measurement and Political
#' Research: An Optimistic View} American Journal of Political Science 43(1):
#' 271-301.
#' 
#' Young, Forrest.  1981.  \sQuote{Quantitative Analysis of Qualitative Data}
#' Psychometrika, 46: 357-388.
#' 
#' Young, Forrest, Jan de Leeuw and Yoshio Takane.  1976.  \sQuote{Regression
#' with Qualitative and Quantitative Variables: An Alternating Least Squares
#' Method with Optimal Scaling Features} Psychometrika, 41:502-529.
balsos <- function(formula, data, iter=2500, chains = 1, 
                   alg = c("NUTS", "HMC", "Fixed_param"), ...){
  if(!requireNamespace("rstan")){
    stop("The rstan package must be installed to use balsos.\n")
  }
  alg <- match.arg(alg)
  stancode <- "data{
  int N; //number of observations
  int M; //number of DV categories
  int k; //number of columns of model matrix
  real ybar; //mean of y
  real sy; //sd of y
  int y[N]; //untransformed DV
  matrix[N,k] X; //model matrix
}
parameters {
  vector[k] beta; //regression coefficients
  real<lower=0> sigma; //standard deviation
  ordered[M] y_star; //transformed dependent variable
}
transformed parameters {
  vector[N] linpred; //linear predictor
  vector[N] y_new; //vector of transformed DV values
  real sd_new; //sd of transformed y
  real mean_new; //mean of transformed y
  vector[N] y_new2; //rescaled y_new;
  linpred = X*beta;
  for(i in 1:N){
    y_new[i] = y_star[y[i]]; //vector of transformed DV values
  }
  sd_new = sd(y_new);
  mean_new = mean(y_new);
  y_new2 = (y_new - mean_new)*(sy/sd_new) + ybar;
}
model{
  beta[1] ~ cauchy(0,10); //prior on intercept
  for(i in 2:k){
    beta[i] ~ cauchy(0, 2.5); //prior on regression coefficients
  }
    target += normal_lpdf(y_new2 | linpred, sigma);
}"

X <- model.matrix(formula, data)
y <- as.integer(model.response(model.frame(formula, data)))
y <- (y - min(y)) + 1L
balsos.dat <- list(
  M=max(y), N=length(y), k = ncol(X), ybar = mean(y), sy = sd(y),
  y=y, X=X)
fit <- rstan::stan(model_code=stancode, data = balsos.dat,
                   iter = iter, chains = chains, algorithm=alg, ...)
ret <- list(fit = fit, y=y, X=X, form=formula)
class(ret) <- "balsos"
return(ret)
}


#' Testing Measurement Level Assumptions.
#' 
#' Uses the algorithm discussed in Armstrong and Jacoby 2018) to test the
#' intervality of a variable scaled by the Bayesian ALSOS method.
#' 
#' 
#' @param obj An object of class \code{balsos}.
#' @param cor.type Type of correlation to be used in the p-value calculation.
#' @returns Printed output giving the Bayesian p-value evaluating the null that
#' ther eis no interesting difference between the original values and the
#' optimally scaled values.
#' 
#' @export
#' 
#' @author Dave Armstrong
test.balsos <- function(obj, cor.type=c("pearson", "spearman")){
  cor.type <- match.arg(cor.type)
  chains <- as.matrix(obj$fit)
  chains <- t(chains[,grep("y_new2", colnames(chains))])
  yvals <- aggregate(1:length(obj$y), list(obj$y), min)
  chains <- chains[yvals[,2], ]
  r1 <- c(cor(chains, yvals[,1], method=cor.type)^2)
  refs <- sample(1:ncol(chains), ncol(chains), replace=FALSE)
  r0 <- diag(cor(chains, chains[,refs])^2)
  cat("Test of Linearity (higher values indicate higher probability of linearity)\n")
  sprintf("%.3f", mean(r1 >= r0))
}

#' Plot Results from BALSOS
#' 
#' Plots the optimally scaled points with posterior 95\% credible intervals.
#' 
#' 
#' @param x Object of class \code{balsos}.
#' @param freq Logical indicating whether you want the frequentist result
#' plotted alongside the Bayesian result.
#' @param offset If \code{freq=T}, the Bayesian points will be plotted at
#' \code{x-offset} and the frequentist points will be plotted at
#' \code{x+offset}.
#' @param ... Other arguments to be passed down, currently not implement and
#' may conflict with the lattice figure.  To change the figure the best advice
#' would be to save the plot as an oject and use the \code{update} function to
#' change its features.
#' @param plot Logical indicating whether the plot should be returned or just the data. 
#' @return A lattice graph produce by a call to \code{xyplot}.
#' @method plot balsos
#' 
#' @export
#' 
#' @author Dave Armstrong
plot.balsos <- function(x, ..., freq=TRUE, offset=.1, plot=TRUE){
  if(freq){
    Xdat <- as.data.frame(x$X[,-1])
    names(Xdat) <- paste("var", 1:ncol(Xdat), sep="")
    Xdat$y <- x$y
    
    a1 <- alsosDV(y ~ ., data=Xdat)
  }
  
  mins <- aggregate(1:length(x$y), list(x$y), min)
  chains <- as.matrix(x$fit)
  chains <- chains[,grep("y_new2", colnames(chains))]
  os <- chains[, mins[,2]]
  
  s<- summary(as.mcmc(os))
  
  newdf <- data.frame(
    freq = a1$result$os[mins[,2]],
    os = s$statistics[,1],
    lower = s$quantiles[,1],
    upper = s$quantiles[,5],
    orig = sort(unique(x$y)), stringsAsFactors=TRUE
  )
  if(!plot){
    return(newdf)
  }else{
    tp <- trellis.par.get()$superpose.symbol
    if(!freq){
      xyplot(os ~ orig, data=newdf,
             panel = function(x,y, ...){
               panel.points(x,y, cex=.6, col=tp$col[1], pch=tp$pch[1])
               panel.segments(x, newdf$lower, x, newdf$upper, col="black", cols=tp$col[1])
             })
    }
    else{
      xyplot(os ~ orig, data=newdf,
             panel = function(x,y, ...){
               panel.points(x-offset,y, cex=.6, col=tp$col[1], pch=tp$pch[1])
               panel.segments(x-offset, newdf$lower, x-offset, newdf$upper, col=tp$col[1])
               panel.points(x+offset, newdf$freq, col=tp$col[2], pch=tp$pch[2])
             })
    }
    
  }
}

##' Summry method for Bayesian ALSOS
##' 
##' @description summary method for objects of class \code{balsos}
##' @param object Object of class \code{balsos}
##' @param ... Other arguments, currently unimplemented
##' 
##' @export
##' 
##' @method summary balsos
summary.balsos <- function(object, ...){
  coef.inds <- grep("^b", names(object$fit))
  mins <- aggregate(1:length(object$y), list(object$y), min)
  os.inds <- sapply(paste("y_new2[", mins[,2], "]", sep=""), function(z)grep(z, names(object$fit), fixed=TRUE))
  names(object$fit)[c(coef.inds, os.inds)] <- c(colnames(object$X),
                                                paste0("y_", sort(unique(object$y))))
  summary(object$fit, pars=c(colnames(object$X), paste0("y_", sort(unique(object$object)))))
}
