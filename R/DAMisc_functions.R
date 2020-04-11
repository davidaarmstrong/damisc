pGumbel <- function (q, mu = 0, sigma = 1){
  stopifnot(sigma > 0)
  exp(-exp(-((q - mu)/sigma)))
}


#' Scalar Measures of Fit for Binary Variable Models
#' 
#' Calculates scalar measures of fit for models with binary dependent variables
#' along the lines described in Long (1997) and Long and Freese (2005).
#' 
#' \code{binfit} calculates scalar measures of fit (many of which are
#' pseudo-R-squared measures) to describe how well a model fits data with a
#' binary dependent variable.
#' 
#' @param mod A model of class \code{glm} with \code{family=binomial}.
#' @return A named vector of scalar measures of fit
#' @author Dave Armstrong
#' @references Long, J.S.  1997.  Regression Models for Categorical and Limited
#' Dependent Variables.  Thousand Oaks, CA: Sage.
#' 
#' Long, J.S. and J. Freese.  2005.  Regression Models for Categorical Outcomes
#' Using Stata, 2nd ed.  College Station, TX: Stata Press.
#' 
#' @export
#' 
#' @examples
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2), data=france, family=binomial)
#' binfit(left.mod)
#' 
binfit <-
function (mod)
{
    mod <- update(mod, x = TRUE, y = TRUE)
    y <- mod$y
    null.mod <- update(mod, ".~1", data=model.frame(mod))
    b <- mod$coef[-1]
    var.ystar <- t(b) %*% var(model.matrix(mod)[, -1]) %*% b
    G <- -2 * (logLik(null.mod) - logLik(mod))
    res.col1 <- c(logLik(null.mod), deviance(mod), NA, 1 - (logLik(mod)/logLik(null.mod)),
        1 - exp(2*(logLik(null.mod) - logLik(mod))/length(mod$residuals)),
        var.ystar/(var.ystar + switch(mod$family[[2]], logit = pi^2/3,
            probit = 1)), mean(mod$y == as.numeric(fitted(mod) >
            0.5)), BIC(mod))
    res.col2 <- c(logLik(mod), G, pchisq(G, 8, lower.tail = FALSE),
        1 - ((logLik(mod) - mod$rank)/logLik(null.mod)), res.col1[5]/(1 -
            (exp(2*logLik(null.mod)/length(mod[["residuals"]])))),
        1 - (sum((y - fitted(mod))^2)/sum((y - mean(y))^2)),
        (sum(mod$y == as.numeric(fitted(mod) > 0.5)) - max(table(y)))/(length(mod$residuals) -
            max(table(y))), AIC(mod))
    res.vec <- c(res.col1, res.col2)[-3]
    res.col1 <- sprintf("%3.3f", res.col1)
    res.col1[3] <- ""
    res.df <- data.frame(Names1 = c("Log-Lik Intercept Only:",
        paste("D(", mod$df.residual, "):", sep = ""), " ", "McFadden's R2:",
        "ML (Cox-Snell) R2:", "McKelvey & Zavoina R2:", "Count R2:",
        "BIC:"), vals1 = res.col1, Names2 = c("Log-Lik Full Model:",
        paste("LR(", mod$rank - 1, "):", sep = ""), "Prob > LR:",
        "McFadden's Adk R2:", "Cragg-Uhler (Nagelkerke) R2:",
        "Efron's R2:", "Adj Count R2:", "AIC:"), vals2 = sprintf("%3.3f",
        res.col2), stringsAsFactors=TRUE)
    return(res.df)
}


#' Generate Spells for Binary Variables
#' 
#' Beck et. al. (1998) identified that binary time-series cross-section data
#' are discrete-time duration data and time dependence can be modeled in a
#' logistic regression by including a flexible function (e.g., cubic spline) of
#' time since the last event as a covariate.  This function creates the
#' variable identifying time since last event.
#' 
#' 
#' @param data A data frame.
#' @param event Character string giving the name of the dichotomous variable
#' identifying the event (where an event is coded 1 and the absence of an event
#' is coded 0).
#' @param tvar Character string giving the name of the time variable.
#' @param csunit Character string giving the name of the cross-sectional unit.
#' @param pad.ts Logical indicating whether the time-series should be filled
#' in, when panels are unbalanced.
#' @return The original data frame with one additional variable.  The
#' \code{spell} variable identifies the number of observed periods since the
#' last event.
#' @author Dave Armstrong
#' @references Alvarez, M., J.A. Cheibub, F. Limongi and A. Przeworski. 1996.
#' Classifying political regimes. Studies in Comparative International
#' Development 31 (Summer): 1-37.
#' 
#' Beck, N.. J. Katz and R. Tucker. 1998. Beyond Ordinary Logit: Taking Time
#' Seriously in Binary-Time-Series-Cross-Section Models. American Journal of
#' Political Science 42(4): 1260-1288.
#' 
#' @export
#' 
#' @examples
#' 
#' library(splines)
#' ## Data from Alvarez et. al. (1996)
#' data(aclp)
#' newdat <- btscs(aclp, "democ", "year", "country")	
#' 
#' # Estimate Model with and without spell
#' full.mod <- glm(democ ~ log(gdpw) + popg + bs(spell, df=4), data=newdat, family=binomial)	
#' restricted.mod <- glm(democ ~ log(gdpw) + popg, data=newdat, family=binomial)	
#' 	
#' # Incremental F-test of time dependence
#' anova(restricted.mod, full.mod, test='Chisq')	
#' 
btscs <-
function (data, event, tvar, csunit, pad.ts = FALSE)
{
    data$orig_order <- 1:nrow(data)
    data <- data[order(data[[csunit]], data[[tvar]]), ]
    spells <- function(x) {
        tmp <- rep(0, length(x))
        runcount <- 0
        for (j in 2:length(x)) {
            if (x[j] == 0 & x[(j - 1)] == 0) {
                tmp[j] <- runcount <- runcount + 1
            }
            if (x[j] != 0 & x[(j - 1)] == 0) {
                tmp[j] <- runcount + 1
                runcount <- 0
            }
            if (x[j] == 0 & x[(j - 1)] != 0) {
                tmp[j] <- runcount <- 0
            }
        }
        tmp
    }
    sp <- split(data, data[[csunit]])
    if (pad.ts) {
        sp <- lapply(sp, function(x) x[match(seq(min(x[[tvar]],
            na.rm = TRUE), max(x[[tvar]], na.rm = TRUE)), x[[tvar]]),
            ])
        for (i in 1:length(sp)) {
            if (any(is.na(sp[[i]][[event]]))) {
                sp[[i]][[event]][which(is.na(sp[[i]][[event]]))] <- 1
            }
            if (any(is.na(sp[[i]][[tvar]]))) {
                sp[[i]][[tvar]] <- seq(min(sp[[i]][[tvar]], na.rm = TRUE),
                  max(sp[[i]][[tvar]], na.rm = TRUE))
            }
            if (any(is.na(sp[[i]][[csunit]]))) {
                sp[[i]][[csunit]][which(is.na(sp[[i]][[csunit]]))] <- mean(sp[[i]][[csunit]],
                  na.rm = TRUE)
            }
        }
    }
    sp <- lapply(1:length(sp), function(x) {
        cbind(sp[[x]], data.frame(spell = spells(sp[[x]][[event]])), stringsAsFactors=TRUE)
    })
    data <- do.call(rbind, sp)
    if (!pad.ts) {
        if (any(is.na(data$orig_order))) {
            data <- data[-which(is.na(data$orig_order)), ]
        }
        data <- data[data$orig_order, ]
    }
    else {
        data <- data[order(data[[csunit]], data[[tvar]]), ]
    }
    invisible(data)
}



#' Test for Combining Categories in Multinomial Logistic Regression Models.
#' 
#' Tests the null hypothesis that categories can be combined in Multinomial
#' Logistic Regression Models
#' 
#' 
#' @param obj An object of class \code{multinom}.
#' @return A matrix of test statistics and p-values.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
#' combTest(mnl.mod)
#' 
#' 
combTest <- function(obj){
	y <- model.response(model.frame(obj))
	b <- c(t(coef(obj)))
	names(b) <- rownames(vcov(obj))
    names(b) <- gsub("[^a-zA-Z0-9\\:.]", "", names(b))
	l <- levels(y)
	l <- gsub("[^a-zA-Z0-9\\:.]", "", l)
	combs <- combn(l, 2)
	res <- array(dim=dim(t(combs)))
	k <- 1
	inds1 <- which(combs[1,] == l[1])
	for(j in 1:length(inds1)){
		inds <- grep(paste("^", combs[2,inds1[j]], ":", sep=""), names(b))
		inds <- inds[-grep("Intercept", names(b)[inds])]
		Q <- matrix(0, ncol=length(b), nrow=length(inds))
		Q[cbind(1:nrow(Q), inds)] <- 1
		w <- t(Q %*% b) %*% solve(Q %*% vcov(obj) %*% t(Q)) %*% Q %*%b
		p <- pchisq(w, nrow(Q), lower.tail=FALSE)
		res[k,]<- c(w,p)
		k <- k+1
	}
	inds2 <- (1:ncol(combs))[-inds1]
	for(j in 1:length(inds2)){
		indsa <- grep(paste("^", combs[2,inds2[j]], ":", sep=""), names(b))
		indsa <- indsa[-grep("Intercept", names(b)[indsa])]
		indsb <- grep(paste("^", combs[1,inds2[j]], ":", sep=""), names(b))
		indsb <- indsb[-grep("Intercept", names(b)[indsb])]
		Q <- matrix(0, ncol=length(b), nrow=length(inds))
		Q[cbind(1:nrow(Q), indsa)] <- 1
		Q[cbind(1:nrow(Q), indsb)] <- -1
		w <- t(Q %*% b) %*% solve(Q %*% vcov(obj) %*% t(Q)) %*% Q %*%b
		p <- pchisq(w, nrow(Q), lower.tail=FALSE)
		res[k,]<- c(w,p)
		k <- k+1
	}
	ns <- max(nchar(as.character(round(res[,1], 0))))
	r1 <- sprintf(paste("%", ns+4, ".3f", sep=""), res[,1])
    r2 <- sprintf("%.3f", res[,2])
    res <- cbind(r1, r2)
	rownames(res) <- apply(combs, 2, paste, collapse="-")
	colnames(res) <- c("Estimate", "p-value")
	print(res, quote=FALSE)
}



#' Surface Plots for Two-Way Interactions
#' 
#' Makes surface plots to display interactions between two continuous variables
#' 
#' This function makes a surface plot of an interaction between two continuous
#' covariates.  If the model is \cr \deqn{y_{i} = b_{0} + b_{1}x_{i1} +
#' b_{2}x_{i2} + b_{3}x_{i1}\times x_{i2} + \ldots + e_{i},}{y = b[0] + b[1]x1
#' + b[2]x2 + b[3]x1*x2 + ... + e[i],}\cr this function plots \eqn{b_{1}x_{i1}
#' + b_{2}x_{i2} + b_{3}x_{i1}\times x_{i2}}{b[1]x1 + b[2]x2 + b[3]x1*x2} for
#' values over the range of \eqn{X_{1}}{X1} and \eqn{X_{2}}{X2}.  The highest
#' 75\%, 50\% and 25\% of the bivariate density of \eqn{X_{1}}{X1} and
#' \eqn{X_{2}}{X2} (as calculated by \code{sm.density} from the \code{sm}
#' package) are colored in with colors of increasing gray-scale.
#' 
#' @param obj A model object of class \code{lm}
#' @param varnames A two-element character vector where each element is the
#' name of a variable involved in a two-way interaction.
#' @param theta Angle defining the azimuthal viewing direction to be passed to
#' \code{persp}
#' @param phi Angle defining the colatitude viewing direction to be passed to
#' \code{persp}
#' @param xlab Optional label to put on the x-axis, otherwise if \code{NULL},
#' it will take the first element of \code{varnames}
#' @param ylab Optional label to put on the y-axis, otherwise if \code{NULL},
#' it will take the second element of \code{varnames}
#' @param zlab Optional label to put on the z-axis, otherwise if \code{NULL},
#' it will be \sQuote{Predictions}
#' @param hcols Vector of four colors to color increasingly high density areas
#' @param ... Other arguments to be passed down to the initial call to
#' \code{persp}
#' @return \item{x1}{Values of the first element of \code{varnames} used to
#' make predictions.} \item{x2}{Values of the second element of \code{varnames}
#' used to make predictions.} \item{pred}{The predictions based on the values
#' \code{x1} and \code{x2}.} \item{graph}{A graph is produced, but no other
#' information is returned.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(InteractionEx)
#' mod <- lm(y ~ x1*x2 + z, data=InteractionEx)
#' DAintfun(mod, c("x1", "x2"))
#' 
DAintfun <-
function (obj, varnames, theta = 45, phi = 10, xlab=NULL, ylab=NULL, zlab=NULL,
    hcols=NULL, ...)
{
    if (length(varnames) != 2) {
        stop("varnames must be a vector of 2 variable names")
    }
	if(!all(varnames %in% names(obj$coef) )){
		stop("not all variables in varnames are in the model")
	}
    else {
        v1 <- varnames[1]
        v2 <- varnames[2]
        ind <- unique(c(grep(v1, names(obj$coef)), grep(v2, names(obj$coef))))
        ind <- ind[order(ind)]
        b <- obj$coef[ind]
        mod.x <- model.matrix(obj)
        not.ind <- c(1:ncol(mod.x))[!(c(1:ncol(mod.x)) %in% ind)]
        mod.x[, not.ind] <- 0
        dens <- kde2d(mod.x[, varnames[1]], mod.x[,varnames[2]])
        b <- obj$coef[ind]
        v1.seq <- dens$x
        v2.seq <- dens$y
        eff.fun <- function(x1, x2) {
            b[1] * x1 + b[2] * x2 + b[3] * x1 * x2
        }
        if(is.null(hcols)){
          hcols <- paste("gray", seq(from = 20, to = 80, length = 4),
            sep = "")
        }
        if(length(hcols) != 4){
          warning("hcols must have 4 values, using gray palette\n")
          hcols <- paste("gray", seq(from = 20, to = 80, length = 4),
            sep = "")
        }
        predsurf <- outer(v1.seq, v2.seq, eff.fun)
        cutoff <- quantile(c(dens$z), prob = c(0.25, 0.5,
            0.75))
        pred1 <- predsurf
        pred1[dens$z < cutoff[1]] <- NA
        pred2 <- predsurf
        pred2[dens$z < cutoff[2]] <- NA
        pred3 <- predsurf
        pred3[dens$z < cutoff[3]] <- NA
        persp(v1.seq, v2.seq, predsurf,
			xlab = ifelse(is.null(xlab), toupper(v1), xlab),
			ylab = ifelse(is.null(ylab), toupper(v2), ylab),
            zlab = ifelse(is.null(zlab), toupper("Predictions"), zlab),
			col = hcols[1], theta = theta, phi = phi,...)
        par(new = TRUE)
            persp(v1.seq, v2.seq, pred1, col = hcols[2], axes = FALSE,
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi,
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq),
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
        par(new = TRUE)
        persp(v1.seq, v2.seq, pred2, col = hcols[3], axes = FALSE,
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi,
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq),
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
        par(new = TRUE)
        persp(v1.seq, v2.seq, pred3, col = hcols[4], axes = FALSE,
            xlab = "", ylab = "", zlab = "", theta = theta, phi = phi,
            zlim = c(min(c(predsurf)), max(c(predsurf))), ylim = c(min(v2.seq),
                max(v2.seq)), xlim = c(min(v1.seq), max(v1.seq)))
	    invisible(list(x1=v1.seq, x2=v2.seq, pred=predsurf))
    }
}


#' Conditional Effects Plots for Interactions in Linear Models
#' 
#' Generates two conditional effects plots for two interacted continuous
#' covariates in linear models.
#' 
#' This function produces graphs along the lines suggested by Brambor, Clark
#' and Golder (2006) and Berry, Golder and Milton (2012), that show the
#' conditional effect of one variable in an interaction given the values of the
#' conditioning variable.  This is an alternative to the methods proposed by
#' John Fox in his \code{effects} package, upon which this function depends
#' heavily.\cr\cr Specifically, if the model is \cr \deqn{y_{i} = b_{0} +
#' b_{1}x_{i1} + b_{2}x_{i2} + b_{3}x_{i1}\times x_{i2} + \ldots + e_{i},}{y =
#' b[0] + b[1]x1 + b[2]x2 + b[3]x1*x2 + ... + e[i],}\cr this function plots
#' calculates the conditional effect of \eqn{X_{1}}{X1} given
#' \eqn{X_{2}}{X2}\cr \deqn{\frac{\partial y}{\partial X_{1}} = b_{1} +
#' b_{3}X_{2}}{dy/dX1 = b[1] + b[3]X2} and the variances of the conditional
#' effects\cr \deqn{V(b_{1} + b_{3}X_{2}) = V(b_{1} + X_{2}^{2}V(b_{3}) +
#' 2(1)(X_{2})V(b_{1},b_{3}))}{V(b[1] + b[3]X[2]) = V(b[1] + (X[2]^2)V(b[3]) +
#' 2(1)(X[2])V(b[1],b[3]))} for different values of \eqn{X_{2}}{X2} and then
#' switches the places of \eqn{X_{1}}{X1} and \eqn{X_{2}}{X2}, calculating the
#' conditional effect of \eqn{X_{2}}{X2} given a range of values of
#' \eqn{X_{1}}{X1}.  95\% confidence bounds are then calculated and plotted for
#' each conditional effects along with a horizontal reference line at 0.
#' 
#' @param obj A model object of class \code{lm}
#' @param varnames A two-element character vector where each element is the
#' name of a variable involved in a two-way interaction.
#' @param varcov A variance-covariance matrix with which to calculate the
#' conditional standard errors.  If \code{NULL}, it is calculated with
#' \code{vcov(obj)}.
#' @param rug Logical indicating whether a rug plot should be included.
#' @param ticksize A scalar indicating the size of ticks in the rug plot (if
#' included) positive values put the rug inside the plotting region and
#' negative values put it outside the plotting region.
#' @param hist Logical indicating whether a histogram of the x-variable should
#' be included in the plotting region.
#' @param hist.col Argument to be passed to \code{polygon} indicating the color
#' of the histogram bins.
#' @param nclass vector of two integers indicating the number of bins in the
#' two histograms, which will be passed to \code{hist}.
#' @param scale.hist A scalar in the range (0,1] indicating how much vertical
#' space in the plotting region the histogram should take up.
#' @param border Argument passed to \code{polygon} indicating how the border of
#' the histogram bins should be printed (\code{NA} for no border).
#' @param name.stem A character string giving filename to which the appropriate
#' extension will be appended
#' @param xlab Optional vector of length two giving the x-labels for the two
#' plots that are generated.  The first element of the vector corresponds to
#' the figure plotting the conditional effect of the first variable in
#' \code{varnames} given the second and the second element of the vector
#' corresponds to the figure plotting the conditional effect of the second
#' variable in \code{varnames} conditional on the first.
#' @param ylab Optional vector of length two giving the y-labels for the two
#' plots that are generated.  The first element of the vector corresponds to
#' the figure plotting the conditional effect of the first variable in
#' \code{varnames} given the second and the second element of the vector
#' corresponds to the figure plotting the conditional effect of the second
#' variable in \code{varnames} conditional on the first.
#' @param plot.type One of \sQuote{pdf}, \sQuote{png}, \sQuote{eps} or
#' \sQuote{screen}, where the one of the first three will produce two graphs
#' starting with \code{name.stem} written to the appropriate file type and the
#' third will produce graphical output on the screen.
#' @return \item{graphs}{Either a single graph is printed on the screen (using
#' \code{par(mfrow=c(1,2))}) or two figures starting with \code{name.stem} are
#' produced where each gives the conditional effect of one variable based on
#' the values of another.}
#' @author Dave Armstrong
#' @references Brambor, T., W.R. Clark and M. Golder.  (2006) Understanding
#' Interaction Models: Improving Empirical Analyses.  Political Analysis 14,
#' 63-82.\cr Berry, W., M. Golder and D. Milton.  (2012) Improving Tests of
#' Theories Positing Interactions.  Journal of Politics.
#' 
#' @export
#' 
#' @examples
#' 
#' data(InteractionEx)
#' mod <- lm(y ~ x1*x2 + z, data=InteractionEx)
#' DAintfun2(mod, c("x1", "x2"), hist=TRUE, scale.hist=.3)
#' 
DAintfun2 <-
function (obj, varnames, varcov=NULL, rug = TRUE, ticksize = -0.03, hist = FALSE,
    hist.col = "gray75", nclass = c(10, 10), scale.hist = 0.5,
    border = NA, name.stem = "cond_eff",
	xlab = NULL, ylab=NULL, plot.type = "screen")
{
    rseq <- function(x) {
        rx <- range(x, na.rm = TRUE)
        seq(rx[1], rx[2], length = 25)
    }
    MM <- model.matrix(obj)
    v1 <- varnames[1]
    v2 <- varnames[2]
    ind1 <- grep(paste0("^",v1,"$"), names(obj$coef))
    ind2 <- grep(paste0("^",v2,"$"), names(obj$coef))
    indboth <- which(names(obj$coef) %in% c(paste0(v1,":",v2),paste0(v2,":",v1)))
    ind1 <- c(ind1, indboth)
    ind2 <- c(ind2, indboth)
    s1 <- rseq(MM[, v1])
    s2 <- rseq(MM[, v2])
    a1 <- a2 <- matrix(0, nrow = 25, ncol = ncol(MM))
    a1[, ind1[1]] <- 1
    a1[, ind1[2]] <- s2
    a2[, ind2[1]] <- 1
    a2[, ind2[2]] <- s1
    eff1 <- a1 %*% obj$coef
    if(is.null(varcov)){
        varcov <- vcov(obj)
    }
    se.eff1 <- sqrt(diag(a1 %*% varcov %*% t(a1)))
    low1 <- eff1 - qt(0.975, obj$df.residual) * se.eff1
    up1 <- eff1 + qt(0.975, obj$df.residual) * se.eff1
    eff2 <- a2 %*% obj$coef
    se.eff2 <- sqrt(diag(a2 %*% varcov %*% t(a2)))
    low2 <- eff2 - qt(0.975, obj$df.residual) * se.eff2
    up2 <- eff2 + qt(0.975, obj$df.residual) * se.eff2
    if (!plot.type %in% c("pdf", "png", "eps", "screen")) {
        print("plot type must be one of - pdf, png or eps")
    }
    else {
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v1, ".pdf", sep = ""),
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v1, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            old.psopts <- ps.options()
            setEPS()
            postscript(paste(name.stem, "_", v1, ".eps", sep = ""))
        }
        if (plot.type == "screen") {
            oldpar <- par()
            par(mfrow = c(1, 2))
        }
        plot(s2, eff1, type = "n", ylim = range(c(low1, up1)),
            xlab = ifelse(is.null(xlab), toupper(v2), xlab[1]), ylab = ifelse(is.null(ylab), paste("Conditional Effect of ",
                toupper(v1), " | ", toupper(v2), sep = ""), ylab[1]))
        if (hist == TRUE) {
            rng <- diff(par()$usr[3:4])
            h2 <- hist(MM[, v2], nclass = nclass[1], plot = FALSE)
            prop2 <- h2$counts/sum(h2$counts)
            plot.prop2 <- (prop2/max(prop2)) * rng * scale.hist +
                par()$usr[3]
            av2 <- pretty(prop2, n = 3)
            axis(4, at = (av2/max(prop2)) * rng * scale.hist +
                par()$usr[3], labels = av2)
            br2 <- h2$breaks
            for (i in 1:(length(br2) - 1)) {
                polygon(x = c(br2[i], br2[(i + 1)], br2[(i +
                  1)], br2[i], br2[i]), y = c(par()$usr[3], par()$usr[3],
                  plot.prop2[i], plot.prop2[i], par()$usr[3]),
                  col = hist.col, border = border)
            }
        }
        if (rug == TRUE) {
            rug(MM[,v2], ticksize = ticksize)
        }
        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        lines(s2, eff1)
        lines(s2, low1, lty = 2)
        lines(s2, up1, lty = 2)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v2, ".pdf", sep = ""),
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v2, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            postscript(paste(name.stem, "_", v2, ".eps", sep = ""))
        }
        plot(s1, eff2, type = "n", ylim = range(c(low2, up2)),
            xlab = ifelse(is.null(xlab), toupper(v1), xlab[2]),
			ylab = ifelse(is.null(ylab), paste("Conditional Effect of ",
                toupper(v2), " | ", toupper(v1), sep = ""), ylab[2]))
        if (hist == TRUE) {
            rng <- diff(par()$usr[3:4])
            h1 <- hist(MM[,v1], nclass = nclass[2], plot = FALSE)
            prop1 <- h1$counts/sum(h1$counts)
            plot.prop1 <- (prop1/max(prop1)) * rng * scale.hist +
                par()$usr[3]
            av1 <- pretty(prop1, n = 3)
            axis(4, at = (av1/max(prop1)) * rng * scale.hist +
                par()$usr[3], labels = av1)
            br1 <- h1$breaks
            for (i in 1:(length(br1) - 1)) {
                polygon(x = c(br1[i], br1[(i + 1)], br1[(i +
                  1)], br1[i], br1[i]), y = c(par()$usr[3], par()$usr[3],
                  plot.prop1[i], plot.prop1[i], par()$usr[3]),
                  col = hist.col, border = border)
            }
        }
        if (rug == TRUE) {
            rug(MM[, v1], ticksize = ticksize)
        }
        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        lines(s1, eff2)
        lines(s1, low2, lty = 2)
        lines(s1, up2, lty = 2)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "eps") {
            ps.options <- old.psopts
        }
        if (plot.type == "screen") {
            par <- oldpar
        }
    }
}


#' Maximal First Differences for Generalized Linear Models
#' 
#' For objects of class \code{glm}, it calculates the change in predicted
#' responses, for maximal discrete changes in all covariates holding all other
#' variables constant at typical values.
#' 
#' The function calculates the changes in predicted responses for maximal
#' discrete changes in the covariates, for objects of class \code{glm}.  This
#' function works with polynomials specified with the \code{poly} function.  It
#' also works with multiplicative interactions of the covariates by virtue of
#' the fact that it holds all other variables at typical values.  By default,
#' typical values are the median for quantitative variables and the mode for
#' factors. The way the function works with factors is a bit different.  The
#' function identifies the two most different levels of the factor and
#' calculates the change in predictions for a change from the level with the
#' smallest prediction to the level with the largest prediction.
#' 
#' @param obj A model object of class \code{glm}.
#' @param data Data frame used to fit \code{object}.
#' @param typical.dat Data frame with a single row containing values at which
#' to hold variables constant when calculating first differences.  These values
#' will be passed to \code{predict}, so factors must take on a single value,
#' but have all possible levels as their levels attribute.
#' @param diffchange A string indicating the difference in predictor values to
#' calculate the discrete change.  \code{range} gives the difference between
#' the minimum and maximum, \code{sd} gives plus and minus one-half standard
#' deviation change around the median and \code{unit} gives a plus and minus
#' one-half unit change around the median.
#' @param sim Logical indicating whether simulated confidence bounds on the
#' difference should be calculated and presented.
#' @param R Number of simulations to perform if \code{sim} is \code{TRUE}
#' @return A list with the following elements: \item{diffs}{A matrix of
#' calculated first differences} \item{minmax}{A matrix of values that were
#' used to calculate the predicted changes}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2), data=france, family=binomial)
#' typical.france <- data.frame(
#' 	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
#' 	age = 35, stringsAsFactors=TRUE
#' 	)
#' glmChange(left.mod, data=france, typical.dat=typical.france)
#' 
glmChange <-
function (obj, data, typical.dat = NULL, diffchange = c("range", "sd", "unit"), sim=FALSE, R=1000)
{
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]),
                sep = "")
            col.inds <- match(tmp.levs, names(obj$coef))
            if (length(grep("1$", names(obj$coef)[col.inds])) >
                0) {
                col.inds <- c(col.inds[which(is.na(col.inds))],
                  col.inds[grep("1$", names(col.inds))])
                names(col.inds) <- gsub("1$", "", names(col.inds))
                col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            }
            tmp.coefs <- obj$coef[col.inds]
            tmp.coefs <- obj$coef[match(tmp.levs, names(obj$coef))]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm],
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)],
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	if(length(vars) > 0){
        mmc <- match.arg(diffchange)
        for (i in 1:length(vars)) {
            if(mmc == "range"){
            minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = TRUE)
            }
            if(mmc == "sd"){
              tmp <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=TRUE)
              tmp[1] <- ifelse(tmp[1] < min(data[[vars[i]]], na.rm=TRUE), min(data[[vars[i]]], na.rm=TRUE), tmp[1])
              tmp[2] <- ifelse(tmp[2] > max(data[[vars[i]]], na.rm=TRUE), max(data[[vars[i]]], na.rm=TRUE), tmp[2])
              minmax[[vars[i]]] <- tmp
            }
            if(mmc == "unit"){
            minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)
            }
            meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE)
        }
    }
    tmp.df <- do.call(data.frame, c(lapply(meds, function(x) rep(x,
        length(meds) * 2)), stringsAsFactors=TRUE))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ",
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1,
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1,
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
    preds <- matrix(predict(obj, newdata = tmp.df, type = "response"),
        ncol = 2, byrow = TRUE)
    diffs <- cbind(preds, apply(preds, 1, diff))
    colnames(diffs) <- c("min", "max", "diff")
    rownames(diffs) <- rn
    minmax.mat <- do.call(data.frame, c(minmax, stringsAsFactors=TRUE))
    minmax.mat <- rbind(do.call(data.frame, c(meds, stringsAsFactors=TRUE)), minmax.mat)
    rownames(minmax.mat) <- c("typical", "min", "max")
	if(sim){
    preds <- predict(obj, newdata = tmp.df, type = "link", se.fit=TRUE)
	res <- sapply(1:R, function(x)apply(matrix(family(obj)$linkinv(with(preds, rnorm(length(preds$fit), fit, se.fit))), ncol=2, byrow=TRUE), 1, diff))
	cis <- t(apply(res, 1, quantile, c(.025, .975)))
	colnames(cis) <- c("lower", "upper")
	diffs <- cbind(diffs, cis)
	}
    ret <- list(diffs = diffs, minmax = minmax.mat)
    class(ret) <- "change"
    return(ret)
}


#' Functions for Estimating Interaction Effects in Logit and Probit Models
#' 
#' Norton and Ai (2003) and Norton, Wang and Ai (2004) discuss methods for
#' calculating the appropriate marginal effects for interactions in binary
#' logit/probit models.  These functions are direct translations of the Norton,
#' Wang and Ai (2004) Stata code.
#' 
#' 
#' @param obj A binary logit or probit model estimated with \code{glm}.
#' @param vars A vector of the two variables involved in the interaction.
#' @param data A data frame used in the call to \code{obj}.
#' @return A list is returned with two elements - \code{byobs} and
#' \code{atment}.  The \code{byobs} result gives the interaction effect
#' evaluated at each observation.  The \code{atmean} element has the marginal
#' effect evaluated at the mean.  Each eleement contains an element \code{int}
#' which is a data frame with the following variable: \item{int_eff}{The
#' correctly calucalted marginal effect.} \item{linear}{The incorrectly
#' calculated marginal effect following the linear model analogy.}
#' \item{phat}{Predicted Pr(Y=1|X).} \item{se_int_eff}{Standard error of
#' \code{int_eff}.} \item{zstat}{The interaction effect divided by its standard
#' error}
#' 
#' The \code{X} element of each returned result is the X-matrix used to
#' generate the result.
#' @author Dave Armstrong
#' @references Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  Computing
#' Interaction Effects and Standard Errors in Logit and Probit Models.  The
#' Stata Journal 4(2): 154-167.\cr
#' 
#' Ai, Chunrong and Edward C. Norton.  2003.  Interaction Terms in Logit and
#' Probit Models.  Economics Letters 80(1): 123-129.
#' 
#' Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  inteff: Computing
#' Interaction Effects and Standard Errors in Logit and Probit Models, Stata
#' Code.
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' mod <- glm(voteleft ~ age*lrself + retnat + male, data=france, family=binomial)
#' out <- intEff(obj=mod, vars=c("age", "lrself"), data=france)
#' out <- out$byobs$int
#' plot(out$phat, out$int_eff, xlab="Predicted Pr(Y=1|X)", 
#' 	ylab = "Interaction Effect")
#' ag <- aggregate(out$linear, list(out$phat), mean)
#' lines(ag[,1], ag[,2], lty=2, col="red", lwd=2)
#' legend("topright", c("Correct Marginal Effect", "Linear Marginal Effect"), 
#' 	pch=c(1, NA), lty=c(NA, 2), col=c("black", "red"), lwd=c(NA, 2), inset=.01)
#' 
intEff <-
function (obj, vars, data)
{
    if (obj$family$family != "binomial")
        stop("intEff only works for binomial GLMs")
    dat.cl <- attr(obj$terms, "dataClasses")[vars]
    if (dat.cl[vars[1]] == "factor" & length(unique(data[[vars[1]]])) ==
        2) {
        vars[1] <- paste(vars[1], colnames(contrasts(data[[vars[1]]])),
            sep = "")
    }
    if (dat.cl[vars[1]] == "factor" & length(unique(data[[vars[1]]])) >
        2) {
        stop("factor variables must only have two unique values, violated by v1")
    }
    if (dat.cl[vars[2]] == "factor" & length(unique(data[[vars[2]]])) ==
        2) {
        vars[2] <- paste(vars[2], colnames(contrasts(data[[vars[2]]])),
            sep = "")
    }
    if (dat.cl[vars[2]] == "factor" & length(unique(data[[vars[2]]])) >
        2) {
        stop("factor variables must only have two unique values, violated by v2")
    }
    v1 <- paste(vars, collapse = ":")
    v2 <- paste(vars[c(2, 1)], collapse = ":")
    inb <- any(c(v1, v2) %in% names(obj$coef))
    X <- model.matrix(obj, obj$model)
    if (!inb)
        stop("Specified variables not interacted in model\n")
    if (v2 %in% names(obj$coef)) {
        vars <- vars[c(2, 1)]
    }
    int.var <- paste(vars, collapse = ":")
    b <- obj$coef
    lens <- apply(X[, vars], 2, function(x) length(unique(x)))
    if (any(dat.cl == "numeric"))
        dat.cl[which(dat.cl == "numeric")] <- "c"
    if (any(dat.cl == "factor"))
        dat.cl[which(dat.cl == "factor")] <- "d"
    if (any(lens == 2))
        dat.cl[which(lens == 2)] <- "d"
    type.int <- paste(sort(dat.cl), collapse = "")
    if (obj$family$link == "logit")
        type.int <- paste("l", type.int, sep = "")
    if (obj$family$link == "probit")
        type.int <- paste("p", type.int, sep = "")
    if (!(obj$family$link %in% c("probit", "logit")))
        stop("Link must be either logit or probit")
    out.dat <- switch(type.int, lcc = logit_cc(obj = obj, int.var = int.var,
        vars = vars, b = b, X = X), lcd = logit_cd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = X), ldd = logit_dd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = X), pcc = probit_cc(obj = obj,
        int.var = int.var, vars = vars, b = b, X = X), pcd = probit_cd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = X), pdd = probit_dd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = X))
    meanX <- matrix(colMeans(X), nrow=1)
    colnames(meanX) <- colnames(X)
    mean.out.dat <- switch(type.int, lcc = logit_cc(obj = obj, int.var = int.var,
        vars = vars, b = b, X = meanX), lcd = logit_cd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = meanX), ldd = logit_dd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = meanX), pcc = probit_cc(obj = obj,
        int.var = int.var, vars = vars, b = b, X = meanX), pcd = probit_cd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = meanX), pdd = probit_dd(obj = obj,
        int.var = int.var, vars = vars, b = b, X = meanX))
    res <- list(byobs = list(int = out.dat, X=X), atmean = list(mean.out.dat, X=meanX))
    invisible(res)
}


#' Functions for Estimating Interaction Effects in Logit and Probit Models
#' 
#' Norton and Ai (2003) and Norton, Wang and Ai (2004) discuss methods for
#' calculating the appropriate marginal effects for interactions in binary
#' logit/probit models.  These functions are direct translations of the Norton,
#' Wang and Ai (2004) Stata code.  These functions are not intended to be
#' called by the user directly, rather they are called as needed by
#' \code{intEff}.
#' 
#' 
#' @aliases logit_cc logit_cd logit_dd probit_cc probit_cd probit_dd
#' @param obj A binary logit or probit model estimated with \code{glm}.
#' @param int.var The name of the interaction variable.
#' @param vars A vector of the two variables involved in the interaction.
#' @param b Coefficients from the \code{glm} object.
#' @param X Model matrix from the \code{glm} object.
#' @return A data frame with the following variable: \item{int_eff}{The
#' correctly calucalted marginal effect.} \item{linear}{The incorrectly
#' calculated marginal effect following the linear model analogy.}
#' \item{phat}{Predicted Pr(Y=1|X).} \item{se_int_eff}{Standard error of
#' \code{int_eff}.} \item{zstat}{The interaction effect divided by its standard
#' error}
#' @author Dave Armstrong
#' @references Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  Computing
#' Interaction Effects and Standard Errors in Logit and Probit Models.  The
#' Stata Journal 4(2): 154-167.\cr
#' 
#' Ai, Chunrong and Edward C. Norton.  2003.  Interaction Terms in Logit and
#' Probit Models.  Economics Letters 80(1): 123-129.
#' 
#' Norton, Edward C., Hua Wang and Chunrong Ai.  2004.  inteff: Computing
#' Interaction Effects and Standard Errors in Logit and Probit Models, Stata
#' Code.
#' 
#' @export
#' 
logit_cc <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- (X %*% b)[,,drop=F]
    phat <- plogis(xb)
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    logitcc <- b[int.var] * phi + (b[vars[1]] + b[int.var] *
        X[, vars[2]]) * (b[vars[2]] + b[int.var] * X[, vars[1]]) *
        phi * (1 - (2 * phat))
    d2f <- phat * (1 - phat) * (1 - 2 * phat)
    d3f <- phat * (1 - phat) * (1 - 6 * phat + 6 * phat^2)
    b1b4x2 <- b[vars[1]] + b[int.var] * X[, vars[2]]
    b2b4x1 <- b[vars[2]] + b[int.var] * X[, vars[1]]
    deriv11 <- b[int.var] * d2f * X[, vars[1]] + b2b4x1 * d2f +
        b1b4x2 * b2b4x1 * d3f * X[, vars[1]]
    deriv22 <- b[int.var] * d2f * X[, vars[2]] + b1b4x2 * d2f +
        b1b4x2 * b2b4x1 * X[, vars[2]] * d3f
    deriv44 <- phi + b[int.var] * d2f * X[, vars[1]] * X[, vars[2]] +
        X[, vars[2]] * b2b4x1 * d2f + X[, vars[1]] * b1b4x2 *
        d2f + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] *
        d3f
    derivcc <- b[int.var] * d2f + b1b4x2 * b2b4x1 * d3f
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) b[int.var] * d2f * x +
        b1b4x2 * b2b4x1 * x * d3f)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv11, deriv22, deriv44, nn, derivcc)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitcc/logit_se
    out <- data.frame(int_eff = logitcc, linear = linear, phat = phat,
        se_int_eff = logit_se, zstat = logit_t, stringsAsFactors=TRUE)
    invisible(out)
}
#' 
#' @export
#' 
logit_cd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- X %*% b
    phat <- plogis(xb)
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    dum <- vars[which(sapply(apply(X[, vars], 2, table), length) ==
        2)]
    cont <- vars[which(vars != dum)]
    X1 <- X2 <- X
    X1[, dum] <- 1
    X1[, int.var] <- X1[, cont] * X1[, dum]
    phat1 <- plogis(X1 %*% b)
    phi1 <- phat1 * (1 - phat1)
    d2f1 <- phi1 * (1 - 2 * phat1)
    ie1 <- (b[cont] + b[int.var]) * phi1
    X2[, dum] <- 0
    X2[, int.var] <- X2[, cont] * X2[, dum]
    phat2 <- plogis(X2 %*% b)
    phi2 <- phat2 * (1 - phat2)
    d2f2 <- phi2 * (1 - 2 * phat2)
    ie2 <- b[cont] * phi2
    logitcd <- ie1 - ie2
    deriv1 <- phi1 - phi2 + b[cont] * X[, cont] * (d2f1 - d2f2) +
        b[int.var] * X[, cont] * d2f1
    deriv2 <- (b[cont] + b[int.var]) * d2f1
    deriv3 <- phi1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
    deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) ((b[cont] + b[int.var]) *
        d2f1 - b[cont] * d2f2) * x)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitcd/logit_se
    out <- data.frame(int_eff = logitcd, linear = linear, phat = phat,
        se_int_eff = logit_se, zstat = logit_t, stringsAsFactors=TRUE)
    invisible(out)
}
#' 
#' @export
#' 
logit_dd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- X %*% b
    phat <- plogis(xb)
    phi <- phat * (1 - phat)
    linear <- b[int.var] * phi
    X11 <- X01 <- X10 <- X00 <- X
    X11[, vars[1]] <- 1
    X11[, vars[2]] <- 1
    X10[, vars[1]] <- 1
    X10[, vars[2]] <- 0
    X01[, vars[1]] <- 0
    X01[, vars[2]] <- 1
    X00[, vars[1]] <- 0
    X00[, vars[2]] <- 0
    X00[, int.var] <- X00[, vars[1]] * X00[, vars[2]]
    X11[, int.var] <- X11[, vars[1]] * X11[, vars[2]]
    X01[, int.var] <- X01[, vars[1]] * X01[, vars[2]]
    X10[, int.var] <- X10[, vars[1]] * X10[, vars[2]]
    phat11 <- plogis(X11 %*% b)
    phat00 <- plogis(X00 %*% b)
    phat10 <- plogis(X10 %*% b)
    phat01 <- plogis(X01 %*% b)
    phi11 <- phat11 * (1 - phat11)
    phi10 <- phat10 * (1 - phat10)
    phi01 <- phat01 * (1 - phat01)
    phi00 <- phat00 * (1 - phat00)
    logitdd <- (phat11 - phat10) - (phat01 - phat00)
    deriv1 <- phi11 - phi10
    deriv2 <- phi11 - phi01
    deriv3 <- phi11
    deriv0 <- (phi11 - phi01) - (phi10 - phi00)
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) ((phi11 - phi01) - (phi10 -
        phi00)) * x)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    logit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    logit_t <- logitdd/logit_se
    out <- data.frame(int_eff = logitdd, linear = linear, phat = phat,
        se_int_eff = logit_se, zstat = logit_t, stringsAsFactors=TRUE)
    invisible(out)
}


#' Print Statistically Significant MNL Coefficients
#' 
#' By default, the summary for objects of class \code{multinom} is not
#' particularly helpful.  It still requires a lot of work on the part of the
#' user to figure out which coefficients are significantly different from zero
#' and which ones are not.  \code{mnlSig} solves this problem by either
#' flagging significant coefficients with an asterisk or only printing
#' significant coefficients, leaving insignificant ones blank.
#' 
#' 
#' @param obj A model object of class \code{multinom}.
#' @param pval The desired Type I error rate to identify coefficients as
#' statistically significant.
#' @param two.sided Logical indicating whether calculated p-values should be
#' two-sided (if \code{TRUE}) or one-sided (if \code{FALSE}).
#' @param flag.sig Logical indicating whether an asterisk should be placed
#' beside coefficients which are significant at the \code{pval} level.
#' @param insig.blank Logical indicating whether coefficients which are not
#' significant at the \code{pval} level should be blank in the output.
#' @return A data frame suitable for printing with the (optionally
#' significance-flagged) coefficients from a multinomial logit model.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ retnat + male + retnat + lrself, data=france)
#' mnlSig(mnl.mod)
#' 
mnlSig <-
function (obj, pval = 0.05, two.sided = TRUE, flag.sig = TRUE,
    insig.blank = FALSE)
{
    smulti <- summary(obj)
    multi.t <- smulti$coefficients/smulti$standard.errors
    multi.p <- (2^as.numeric(two.sided)) * pnorm(abs(multi.t),
        lower.tail = FALSE)
    b <- matrix(sprintf("%.3f", smulti$coefficients), ncol = ncol(multi.t))
    sig.vec <- c(" ", "*")
    sig.obs <- as.numeric(multi.p < pval) + 1
    if (flag.sig) {
        b <- matrix(paste(b, sig.vec[sig.obs], sep = ""), ncol = ncol(multi.t))
    }
    if (insig.blank) {
        b[which(multi.p > pval, arr.ind = TRUE)] <- ""
    }
    rownames(b) <- rownames(multi.t)
    colnames(b) <- colnames(multi.t)
    b <- as.data.frame(b)
    return(b)
}



#' Average Effects Plot for Multinomial Logistic Regression
#' 
#' Produces a plot of average effects for one variable while holding the others
#' constant at observed values.
#' 
#' 
#' @param obj An object of class \code{multinom}.
#' @param varname A string indicating the variable for which the plot is
#' desired.
#' @param data The data used to estimate \code{obj}.
#' @param R Number of simulations used to generate confidence bounds.
#' @param nvals Number of evaluation points for the predicted probabilities.
#' @param plot Logical indicating whether a plot should be produced (if
#' \code{TRUE}) or numerical results should be returned (if \code{FALSE}).
#' @param \dots Other arguments to be passed down to \code{xyplot}.
#' @return Either a plot or a data frame with variables \item{mean}{The average
#' effect (i.e., predicted probability)} \item{lower}{The lower 95\% confidence
#' bound} \item{upper}{The upper 95\% confidence bound} \item{y}{The values of
#' the dependent variable being predicted} \item{x}{The values of the
#' independent variable being manipulated}
#' @author Dave Armstrong
#' @references Hanmer, M.J. and K.O. Kalkan.  2013.  \sQuote{Behind the Curve:
#' Clarifying the Best Approach to Calculating Predicted Probabilities and
#' Marginal Effects from Limited Dependent Variable Models}.  American Journal
#' of Political Science.  57(1): 263-277.
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
#' \dontrun{mnlAveEffPlot(mnl.mod, "lrself", data=france)}
#' 
mnlAveEffPlot <- function(obj, varname, data, R=1500, nvals=25, plot=TRUE,...){
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
    d0 <- list()
    if(is.numeric(data[[varname]])){
	  s <- seq(min(data[[varname]], na.rm=TRUE), max(data[[varname]], na.rm=TRUE), length=nvals)
	  for(i in 1:length(s)){
		d0[[i]] <- data
		d0[[i]][[varname]] <- s[i]
	   }
	}
    if(!is.numeric(data[[varname]])){
     s <- obj$xlevels[[varname]]
      for(j in 1:length(s)){
        d0[[j]] <- data
        d0[[j]][[varname]] <- factor(j, levels=1:length(s), labels=s)
      }
    }
	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x))
	y <- model.response(model.frame(obj))
	ylev <- levels(y)
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
	xb <- lapply(Xmats, function(x)lapply(1:nrow(b), function(z)cbind(1, exp(x %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=TRUE))))))
	probs <- lapply(xb, function(x)lapply(x, function(z)z/rowSums(z)))
	out.ci <- lapply(probs, function(x)sapply(x, colMeans))
	out.ci <- lapply(out.ci, apply, 1, quantile, c(.5,.025,.975))
	tmp <- data.frame(
		mean = do.call("c", lapply(out.ci, function(x)x[1,])),
		lower = do.call("c", lapply(out.ci, function(x)x[2,])),
		upper = do.call("c", lapply(out.ci, function(x)x[3,])),
		y = rep(ylev, length(out.ci)), stringsAsFactors=TRUE
		)
	if(is.numeric(data[[varname]])){tmp$s <- rep(s, each=length(ylev))}
	else{tmp$s <- factor(rep(1:length(s), each=length(ylev)), labels=s)}
		if(plot){
			if(is.factor(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				scales=list(x=list(at=1:length(s), labels=s)),
				panel = function(x,y, lower, upper, subscripts){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower[subscripts], x, upper[subscripts], lty=1, col="black")
			})
			}
			if(is.numeric(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				panel=panel.ci, zl=FALSE)
			}
			return(pl)
		}
		else{
			return(tmp)
		}
	}




#' Maximal First Differences for Multinomial Logistic Regression Models
#' 
#' For objects of class \code{multinom}, it calculates the change in predicted
#' probabilities, for maximal discrete changes in all covariates holding all
#' other variables constant at typical values.
#' 
#' The function calculates the changes in predicted probabilities for maximal
#' discrete changes in the covariates for objects of class \code{multinom}.
#' This function works with polynomials specified with the \code{poly}
#' function.  It also works with multiplicative interactions of the covariates
#' by virtue of the fact that it holds all other variables at typical values.
#' By default, typical values are the median for quantitative variables and the
#' mode for factors. The way the function works with factors is a bit
#' different.  The function identifies the two most different levels of the
#' factor and calculates the change in predictions for a change from the level
#' with the smallest prediction to the level with the largest prediction.
#' 
#' @param obj A model object of class \code{multinom}.
#' @param data Data frame used to fit \code{object}.
#' @param typical.dat Data frame with a single row containing values at which
#' to hold variables constant when calculating first differences.  These values
#' will be passed to \code{predict}, so factors must take on a single value,
#' but have all possible levels as their levels attribute.
#' @param diffchange A string indicating the difference in predictor values to
#' calculate the discrete change.  \code{range} gives the difference between
#' the minimum and maximum, \code{sd} gives plus and minus one-half standard
#' deviation change around the median and \code{unit} gives a plus and minus
#' one-half unit change around the median.
#' @param sim Logical indicating whether simulated confidence bounds should be
#' produced.
#' @param R Number of simulations to perform if \code{sim = TRUE}
#' @return A list with the following elements: \item{diffs}{A matrix of
#' calculated first differences} \item{minmax}{A matrix of values that were
#' used to calculate the predicted changes} \item{minPred}{A matrix of
#' predicted probabilities when each variable is held at its minimum value, in
#' turn.} \item{maxPred}{A matrix of predicted probabilities when each variable
#' is held at its maximum value, in turn.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
#' typical.france <- data.frame(
#' 	age = 35, 
#' 	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
#' 	stringsAsFactors=TRUE)
#' mnlChange(mnl.mod, data=france, typical.dat=typical.france)	
#' 
mnlChange <-
function (obj, data, typical.dat = NULL, diffchange=c("range", "sd", "unit"),
 	sim=TRUE, R=1500){
	y <- model.response(model.frame(obj))
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]),
                sep = "")
            col.inds <- match(tmp.levs, obj$coef)
            # if (length(grep("1$", obj$coef[col.inds])) >
            #     0) {
            #     col.inds <- c(col.inds[which(is.na(col.inds))],
            #       col.inds[grep("1$", names(col.inds))])
            #     names(col.inds) <- gsub("1$", "", names(col.inds))
            #     col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            # }
            tmp.coefs <- coef(obj)[,col.inds]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(colMeans(tmp.coefs)), which.max(colMeans(tmp.coefs)))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm],
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)],
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	mmc <- match.arg(diffchange)
    for (i in 1:length(vars)) {
		if(mmc == "range"){
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = TRUE)
		}
		if(mmc == "sd"){
		  tmp <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=TRUE)
		  tmp[1] <- ifelse(tmp[1] < min(data[[vars[i]]], na.rm=TRUE), min(data[[vars[i]]], na.rm=TRUE), tmp[1])
		  tmp[2] <- ifelse(tmp[2] > max(data[[vars[i]]], na.rm=TRUE), max(data[[vars[i]]], na.rm=TRUE), tmp[2])
		  minmax[[vars[i]]] <- tmp
		}
		if(mmc == "unit"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)
		}
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE)
    }
    tmp.df <- do.call(data.frame, c(lapply(meds, function(x) rep(x,
        length(meds) * 2)), stringsAsFactors=TRUE))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ",
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1,
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1,
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
	if(!sim){
	    preds <- predict(obj, newdata = tmp.df, type = "probs")
	    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
	    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
	    diffs <- preds.max - preds.min
	    rownames(preds.min) <- rownames(preds.max) <- rownames(diffs) <- rn
	}
	if(sim){
    preds <- predict(obj, newdata = tmp.df, type = "probs")
    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
	X <- model.matrix(as.formula(paste("~", as.character(formula(obj))[3], sep="")), data=tmp.df)
	xb <- lapply(1:nrow(b), function(z)cbind(1, exp(X %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=TRUE)))))
	probs <- lapply(xb, function(x)t(apply(x, 1, function(z)z/sum(z))))
	pminlist <- lapply(probs, function(x)x[seq(1, nrow(probs[[1]]), by=2), ])
	pmaxlist <- lapply(probs, function(x)x[seq(2, nrow(probs[[1]]), by=2), ])
	difflist <- lapply(1:length(probs), function(i)pmaxlist[[i]]-pminlist[[i]])
	eg <- as.matrix(expand.grid(1:nrow(difflist[[1]]), 1:ncol(difflist[[1]])))
	means <- matrix(sapply(1:nrow(eg), function(ind)mean(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]))), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))
	lower <- matrix(sapply(1:nrow(eg), function(ind)quantile(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]), .025)), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))
	upper <- matrix(sapply(1:nrow(eg), function(ind)quantile(sapply(difflist, function(x)x[eg[ind,1], eg[ind,2]]), .975)), ncol=ncol(difflist[[1]]), nrow=nrow(difflist[[1]]))

	colnames(means) <- colnames(lower) <- colnames(upper) <- colnames(preds.min) <- colnames(preds.max) <- levels(y)
	rownames(means) <- rownames(lower) <- rownames(upper) <- rownames(preds.min) <- rownames(preds.max) <- colnames(tmp.df)
	diffs <- list(mean = means, lower=lower, upper=upper)
}
minmax.mat <- do.call(data.frame, c(minmax, stringsAsFactors=TRUE))
minmax.mat <- rbind(do.call(data.frame, c(meds, stringsAsFactors=TRUE)), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")

ret <- list(diffs = diffs, minmax = minmax.mat, minPred = preds.min,
    maxPred = preds.max)
class(ret) <- "ordChange"
return(ret)
}



#' Average Effects for Multinomial Logistic Regression Models
#' 
#' Calculates average effects of a variable in multinomial logistic regression
#' holding all other variables at observed values.
#' 
#' 
#' @param obj An object of class \code{multinom}
#' @param varnames A string identifying the variable to be manipulated.
#' @param data Data frame used to fit \code{object}.
#' @param diffchange A string indicating the difference in predictor values to
#' calculate the discrete change. \code{sd} gives plus and minus one-half
#' standard deviation change around the median and \code{unit} gives a plus and
#' minus one-half unit change around the median.
#' @param R Number of simulations.
#' @return A list with elements: \item{mean}{Average effect of the variable for
#' each category of the dependent variable.} \item{lower}{Lower 95 percent
#' confidence bound} \item{upper}{Upper 95 percent confidence bound}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
#' mnlChange2(mnl.mod, "lrself", data=france, )	
#' 
#' 
#' 
mnlChange2 <-
  function (obj, varnames, data, diffchange=c("unit", "sd"), R=1500)
  {
      vars <- all.vars(formula(obj))[-1]
      if(any(!(vars %in% names(data)))){
          vars <- vars[-which(!vars %in% names(data))]
      }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(t(coef(obj))), vcov(obj))
    change <- match.arg(diffchange)
    allmean <- alllower <- allupper <- NULL
    for(m in 1:length(varnames)){
        if(!is.factor(data[[varnames[m]]])){
          delt <- switch(change,
                         unit=1,
                         sd = sd(data[[varnames[m]]], na.rm=TRUE))
        } else{
          delt <- 1
        }
        d0 <- list()
        if(is.numeric(data[[varnames[m]]])){
          d0[[1]] <- d0[[2]] <- data
          tmp0 <- d0[[1]][[varnames[m]]]-(.5*delt)
          tmp1 <- d0[[2]][[varnames[m]]]+(.5*delt)
          tmp0 <- ifelse(tmp0 < min(d0[[1]][[varnames[m]]], na.rm=TRUE), 
                         min(d0[[1]][[varnames[m]]], na.rm=TRUE), tmp0)
          tmp1 <- ifelse(tmp1 > max(d0[[2]][[varnames[m]]], na.rm=TRUE), 
                         max(d0[[2]][[varnames[m]]], na.rm=TRUE), tmp1)
          d0[[1]][[varnames[m]]] <- tmp0
          d0[[2]][[varnames[m]]] <- tmp1
        }
        if(!is.numeric(data[[varnames[m]]])){
          l <- obj$xlevels[[varnames[m]]]
          for(j in 1:length(l)){
            d0[[j]] <- data
            d0[[j]][[varnames[m]]] <- factor(j, levels=1:length(l), labels=l)
          }
        }
    	y <- model.response(model.frame(obj))
    	ylev <- levels(y)
    	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x))
    	xb <- lapply(Xmats, function(x)lapply(1:nrow(b), function(z)cbind(1, exp(x %*% t(matrix(c(t(b[z,])), ncol=ncol(coef(obj)), byrow=TRUE))))))
    	probs <- lapply(xb, function(x)lapply(x, function(z)z/rowSums(z)))

    	if(is.numeric(data[[varnames[m]]])){
    	diffs <- lapply(1:R, function(x)probs[[2]][[x]] - probs[[1]][[x]])
    	probdiffs <- sapply(diffs, colMeans)

    	pwdiffmean <- apply(probdiffs, 1, quantile, c(.5,.025,.975))
    	means <- matrix(pwdiffmean[1,], ncol=1)
    	lower <- matrix(pwdiffmean[2,], ncol=1)
    	upper <- matrix(pwdiffmean[3,], ncol=1)

    	}
    	if(is.factor(data[[varnames[m]]])){
    	combs <- combn(length(probs), 2)
    	pwdiffprob <- list()
    	for(j in 1:ncol(combs)){
    		pwdiffprob[[j]] <- lapply(1:R, function(i)probs[[combs[2,j]]][[i]] - probs[[combs[1,j]]][[i]])
    	}

    	out <- lapply(pwdiffprob, function(x)sapply(x, colMeans))
    	out.ci <- lapply(out, apply, 1, quantile, c(.5,.025,.975))
    	means <- sapply(out.ci, function(x)x[1,])
    	lower <- sapply(out.ci, function(x)x[2,])
    	upper <- sapply(out.ci, function(x)x[3,])
    	}
    	l <- obj$xlevels[[varnames[m]]]
    	if(is.numeric(data[[varnames[m]]])){cn <- varnames[m]}
    	else{
    		cn <- paste(varnames[m], ": ", l[combs[2,]], "-", l[combs[1,]], sep="")
    	}
    	colnames(means) <- colnames(lower) <- colnames(upper) <- cn
    	rownames(means) <- rownames(lower) <- rownames(upper) <- ylev
        allmean <- cbind(allmean, means)
        alllower <- cbind(alllower, lower)
        allupper <- cbind(allupper, upper)

    }
	res <- list(diffs=list(mean = t(allmean), lower=t(alllower), upper=t(allupper)))
    class(res) <- "ordChange"
    res
}




#' Maximal First Differences for Proportional Odds Logistic Regression Models
#' 
#' For objects of class \code{polr}, it calculates the change in predicted
#' probabilities, for maximal discrete changes in all covariates holding all
#' other variables constant at typical values.
#' 
#' The function calculates the changes in predicted probabilities for maximal
#' discrete changes in the covariates for objects of class \code{polr}.  This
#' function works with polynomials specified with the \code{poly} function.  It
#' also works with multiplicative interactions of the covariates by virtue of
#' the fact that it holds all other variables at typical values.  By default,
#' typical values are the median for quantitative variables and the mode for
#' factors. The way the function works with factors is a bit different.  The
#' function identifies the two most different levels of the factor and
#' calculates the change in predictions for a change from the level with the
#' smallest prediction to the level with the largest prediction.
#' 
#' @aliases ordChange 
#' 
#' @param obj A model object of class \code{polr}.
#' @param data Data frame used to fit \code{object}.
#' @param typical.dat Data frame with a single row containing values at which
#' to hold variables constant when calculating first differences.  These values
#' will be passed to \code{predict}, so factors must take on a single value,
#' but have all possible levels as their levels attribute.
#' @param diffchange A string indicating the difference in predictor values to
#' calculate the discrete change.  \code{range} gives the difference between
#' the minimum and maximum, \code{sd} gives plus and minus one-half standard
#' deviation change around the median and \code{unit} gives a plus and minus
#' one-half unit change around the median.
#' @param sim Logical indicating whether or not simulations should be done to
#' generate confidence intervals for the difference.
#' @param R Number of simulations.
#' @return A list with the following elements: \item{diffs}{A matrix of
#' calculated first differences} \item{minmax}{A matrix of values that were
#' used to calculate the predicted changes} \item{minPred}{A matrix of
#' predicted probabilities when each variable is held at its minimum value, in
#' turn.} \item{maxPred}{A matrix of predicted probabilities when each variable
#' is held at its maximum value, in turn.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(MASS)
#' data(france)
#' polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
#' typical.france <- data.frame(
#' 	age = 35, 
#' 	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
#' 	stringsAsFactors=TRUE)
#' ordChange(polr.mod, data=france, typical.dat=typical.france, sim=FALSE)	
#' 
ordChange <-
function (obj, data, typical.dat = NULL, diffchange=c("range", "sd", "unit"),
 	sim=TRUE, R=1500){
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
	probfun <- function(x){
		c(cbind(matrix(x, ncol=(ncol(dmat)-1)), 1) %*% dmat)
	}
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$xlevels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]),
                sep = "")
            col.inds <- match(tmp.levs, names(obj$coef))
            # if (length(grep("1$", names(obj$coef)[col.inds])) >
            #     0) {
            #     col.inds <- c(col.inds[which(is.na(col.inds))],
            #       col.inds[grep("1$", names(obj$coef))])
            #     names(col.inds) <- gsub("1$", "", names(col.inds))
            #     col.inds <- col.inds[match(tmp.levs, names(col.inds))]
            # }
            tmp.coefs <- obj$coef[col.inds]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm],
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)],
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
	mmc <- match.arg(diffchange)
    for (i in 1:length(vars)) {
		if(mmc == "range"){
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = TRUE)
		}
		if(mmc == "sd"){
		    tmp <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)*sd(data[[vars[i]]], na.rm=TRUE)
		    tmp[1] <- ifelse(tmp[1] < min(data[[vars[i]]], na.rm=TRUE), min(data[[vars[i]]], na.rm=TRUE), tmp[1])
		    tmp[2] <- ifelse(tmp[2] > max(data[[vars[i]]], na.rm=TRUE), max(data[[vars[i]]], na.rm=TRUE), tmp[2])
		    minmax[[vars[i]]] <- tmp
		}
		if(mmc == "unit"){
        minmax[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE) + c(-.5,.5)
		}
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE)
    }
    tmp.df <- do.call(data.frame, c(lapply(meds, function(x) rep(x,
        length(meds) * 2)), stringsAsFactors=TRUE))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ",
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1,
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1,
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
	if(!sim){
	    preds <- predict(obj, newdata = tmp.df, type = "probs")
	    preds.min <- as.matrix(preds[seq(1, nrow(preds), by = 2), ])
	    preds.max <- as.matrix(preds[seq(2, nrow(preds), by = 2), ])
	    diffs <- preds.max - preds.min
	    rownames(preds.min) <- rownames(preds.max) <- rownames(diffs) <- rn
	}
	if(sim){
        if("polr" %in% class(obj)){
            b <- mvrnorm(R, c(-coef(obj), obj$zeta), vcov(obj))
        }
        if("clm" %in% class(obj)){
            zeta.ind <- grep("|", names(coef(obj)), fixed=TRUE)
            b.ind <- (1:length(coef(obj)))[-zeta.ind]
            b <- mvrnorm(R, c(-coef(obj)[b.ind], coef(obj)[zeta.ind]),
                vcov(obj)[c(b.ind, zeta.ind),c(b.ind, zeta.ind)])
        }
        if(!(class(obj) %in% c("clm", "polr"))){stop("Simulation requires either a clm or polr object\n")}

    X <- model.matrix(update(formula(obj), NULL ~ .), data=tmp.df)
	intlist <- list()
    if(class(obj) == "polr"){
        ylev <- obj$lev
    }
    if(class(obj) == "clm"){
        ylev <- obj$y.levels
    }

    	for(i in 1:(length(ylev)-1)){
			intlist[[i]] <- matrix(0, ncol=(length(ylev)-1), nrow=nrow(X))
			intlist[[i]][,i] <- 1
		}
	X <- X[,-1]
	tmp <- do.call(rbind, lapply(intlist, function(y)cbind(X,y)))
	cprobs <- plogis(tmp %*% t(b))

	dmat <- matrix(0, ncol=length(ylev), nrow=length(ylev))
	dmat[1,1] <- 1
	for(j in 2:length(ylev)){
		dmat[(j-1), j] <- -1
		dmat[j,j] <- 1
	}
	probs <- t(apply(cprobs, 2, probfun))
	cats <- rep(1:length(ylev), each=nrow(tmp.df))
	problist <- lapply(1:max(cats), function(x)probs[, which(cats == x)])
	pminlist <- lapply(problist, function(x)x[,seq(1, ncol(problist[[1]]), by=2)])
	pmaxlist <- lapply(problist, function(x)x[,seq(2, ncol(problist[[1]]), by=2)])
	difflist <- lapply(1:length(problist), function(i)pmaxlist[[i]]-pminlist[[i]])
	means <- sapply(difflist, colMeans)
	lower <- sapply(difflist, apply, 2, quantile, .025)
	upper <- sapply(difflist, apply, 2, quantile, .975)
	rownames(means) <- rownames(lower) <- rownames(upper) <- colnames(tmp.df)
	colnames(means) <- colnames(lower) <- colnames(upper) <- ylev
	preds.min <- sapply(pminlist, colMeans)
	preds.max <- sapply(pmaxlist, colMeans)
    rownames(preds.min) <- rownames(preds.max)  <- colnames(tmp.df)
	colnames(preds.min) <- colnames(preds.max) <- ylev
	diffs <- list(mean = means, lower=lower, upper=upper)
}
minmax.mat <- do.call(data.frame, c(minmax, stringsAsFactors=TRUE))
minmax.mat <- rbind(do.call(data.frame, c(meds, stringsAsFactors=TRUE)), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")

ret <- list(diffs = diffs, minmax = minmax.mat, minPred = preds.min,
    maxPred = preds.max)
class(ret) <- "ordChange"
print(ret)
invisible(ret)
}



#' Average Effects for Proportional Odds Logistic Regression Models
#' 
#' For objects of class \code{polr}, it calculates the average change in
#' predicted probabilities, for discrete changes in a covariate holding all
#' other variables at their observed values.
#' 
#' The function calculates the changes in predicted probabilities for maximal
#' discrete changes in the covariates for objects of class \code{polr}.  This
#' function works with polynomials specified with the \code{poly} function.  It
#' also works with multiplicative interactions of the covariates by virtue of
#' the fact that it holds all other variables at typical values.  By default,
#' typical values are the median for quantitative variables and the mode for
#' factors. The way the function works with factors is a bit different.  The
#' function identifies the two most different levels of the factor and
#' calculates the change in predictions for a change from the level with the
#' smallest prediction to the level with the largest prediction.
#' 
#' @param obj A model object of class \code{polr}.
#' @param varnames A vector of strings identifying the variable to be
#' manipulated.
#' @param data Data frame used to fit \code{object}.
#' @param diffchange A string indicating the difference in predictor values to
#' calculate the discrete change. \code{sd} gives plus and minus one-half
#' standard deviation change around the median and \code{unit} gives a plus and
#' minus one-half unit change around the median.
#' @param R Number of simulations.
#' @return A list with the following elements: \item{diffs}{A matrix of
#' calculated first differences} \item{minmax}{A matrix of values that were
#' used to calculate the predicted changes} \item{minPred}{A matrix of
#' predicted probabilities when each variable is held at its minimum value, in
#' turn.} \item{maxPred}{A matrix of predicted probabilities when each variable
#' is held at its maximum value, in turn.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(MASS)
#' data(france)
#' polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
#' typical.france <- data.frame(
#' 	age = 35, 
#' 	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
#' 	stringsAsFactors=TRUE)
#' ordChange2(polr.mod, "age", data=france, diffchange="sd")	
#' 
ordChange2 <- function (obj, varnames, data, diffchange=c("sd", "unit"),
      R=1500){
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    if("polr" %in% class(obj)){
        b <- mvrnorm(R, c(-coef(obj), obj$zeta), vcov(obj))
    }
    if("clm" %in% class(obj)){
        zeta.ind <- grep("|", names(coef(obj)), fixed=TRUE)
        b.ind <- (1:length(coef(obj)))[-zeta.ind]
        b <- mvrnorm(R, c(-coef(obj)[b.ind], coef(obj)[zeta.ind]),
            vcov(obj)[c(b.ind, zeta.ind),c(b.ind, zeta.ind)])
    }
    if(!(class(obj) %in% c("clm", "polr"))){stop("Simulation requires either a clm or polr object\n")}
    change <- match.arg(diffchange)
    allmean <- alllower <- allupper <- NULL

    for(m in 1:length(varnames)){
      if(!is.factor(data[[varnames[m]]])){
        delt <- switch(change,
                       unit=1,
                       sd = sd(data[[varnames[m]]], na.rm=TRUE))
      } else{
        delt <- 1
      }
        d0 <- list()
        if(is.numeric(data[[varnames[m]]])){
          d0[[1]] <- d0[[2]] <- data
          tmp0 <- d0[[1]][[varnames[m]]]-(.5*delt)
          tmp1 <- d0[[2]][[varnames[m]]]+(.5*delt)
          tmp0 <- ifelse(tmp0 < min(d0[[1]][[varnames[m]]], na.rm=TRUE), 
                         min(d0[[1]][[varnames[m]]], na.rm=TRUE), tmp0)
          tmp1 <- ifelse(tmp1 > max(d0[[2]][[varnames[m]]], na.rm=TRUE), 
                         max(d0[[2]][[varnames[m]]], na.rm=TRUE), tmp1)
          d0[[1]][[varnames[m]]] <- tmp0
          d0[[2]][[varnames[m]]] <- tmp1
        }
        if(!is.numeric(data[[varnames[m]]])){
          l <- obj$xlevels[[varnames[m]]]
          for(j in 1:length(l)){
            d0[[j]] <- data
            d0[[j]][[varnames[m]]] <- factor(j, levels=1:length(l), labels=l)
          }
        }
    	Xmats <- lapply(d0, function(x)model.matrix(formula(obj), data=x)[,-1])
    	intlist <- list()
            if("polr" %in% class(obj)){
                ylev <- obj$lev
            }
            if("clm" %in% class(obj)){
                ylev <- obj$y.levels
            }
        		for(i in 1:(length(ylev)-1)){
    			intlist[[i]] <- matrix(0, ncol=(length(ylev)-1), nrow=nrow(Xmats[[1]]))
    			intlist[[i]][,i] <- 1
    		}
    	tmp <- lapply(Xmats, function(x)do.call(rbind, lapply(intlist, function(y)cbind(x,y))))
    	cprobs <- lapply(tmp, function(x)cbind(plogis(x %*% t(b))))

    	dmat <- matrix(0, ncol=length(ylev), nrow=length(ylev))
    	dmat[1,1] <- 1
    	for(j in 2:length(ylev)){
    		dmat[(j-1), j] <- -1
    		dmat[j,j] <- 1
    	}
    	probfun <- function(x){
    		c(cbind(matrix(x, ncol=(ncol(dmat)-1)), 1) %*% dmat)
    	}
    	probs <- lapply(cprobs, apply, 2, probfun)
    	combs <- combn(length(probs), 2)
    	pwdiffprob <- list()
    	for(j in 1:ncol(combs)){
    		pwdiffprob[[j]] <- probs[[combs[2,j]]] - probs[[combs[1,j]]]
    	}
    	pwdiffmean <- sapply(pwdiffprob, rowMeans)
    	means <- apply(pwdiffmean, 2, function(x)colMeans(matrix(x, ncol=length(ylev))))
    	lower <- apply(pwdiffmean, 2, function(x)apply(matrix(x, ncol=length(ylev)), 2, quantile, .025))
     	upper <- apply(pwdiffmean, 2, function(x)apply(matrix(x, ncol=length(ylev)), 2, quantile, .975))
    	if(is.numeric(data[[varnames[m]]])){cn <- varnames[m]}
    	else{
    		cn <- paste(varnames[m], ": ", l[combs[2,]], "-", l[combs[1,]], sep="")
    	}
    	colnames(means) <- colnames(lower) <- colnames(upper) <- cn
    	rownames(means) <- rownames(lower) <- rownames(upper) <- ylev
        allmean <- cbind(allmean, means)
        alllower <- cbind(alllower, lower)
        allupper <- cbind(allupper, upper)
    }
	res <- list(diffs=list(mean = t(allmean), lower=t(alllower), upper=t(allupper)))
    class(res) <- "ordChange"
    res
}

##' Print method for ordChange objects
##' 
##' @description Print methods for objects of class \code{ordChange}
##' @param x Object of class \code{ordChange}
##' @param ... Other arguments to be passed down to the function
##' @param digits Number of digits to print
##' @export
##' @method print ordChange
print.ordChange <- function(x, ..., digits=3){
    diffs <- x$diffs
    if("list" %in% class(diffs)){
        sig <- ifelse(sign(diffs$lower) == sign(diffs$upper), "*", " ")
        leads <- ifelse(sign(diffs$mean) == -1, "", " ")
        out <- array(paste(leads, sprintf(paste("%0.", digits, "f", sep=""), diffs$mean), sig, sep=""), dim=dim(diffs$mean))
        dimnames(out) <- dimnames(diffs$mean)
    }
    else{
        out <- array(sprintf(paste("%0.", digits, "f", sep=""), diffs), dim=dim(diffs))
        dimnames(out) <- dimnames(diffs)
    }
    print(out, quote=FALSE)
}




#' Fit Statistics and Specification Test for Multinomial Logistic Regression
#' 
#' Provides fit statistics (pseudo R-squared values) and the Fagerland, Hosmer
#' and Bonfi (2008) specification test for Multinomial Logistic Regression
#' models.
#' 
#' 
#' @aliases mnlfit print.mnlfit
#' @param obj An object of class \code{multinom}
#' @param permute Logical indicating whether to check all base categories for
#' the Fagerland et. al. specification test.
#' @return A list with elements: \item{result}{Fit statistics.}
#' \item{permres}{The results of the base category permutation exercise.}
#' @author Dave Armstrong
#' @references Fagerland M. W., D. W. Hosmer and A. M. Bonfi.  2008.
#' \sQuote{Multinomial goodness-of-fit tests for logistic regression models.}
#' Statistics in Medicine.  27(21): 4238-4253.
#' 
#' @export
#' 
#' @examples
#' 
#' library(nnet)
#' data(france)
#' mnl.mod <- multinom(vote ~ age + male + retnat + lrself, data=france)
#' mnlfit(mnl.mod)
#' 
mnlfit <- function(obj, permute=FALSE){
	obj <- update(obj, trace=FALSE)
	y <- model.response(model.frame(obj))
	pp <- predict(obj, type="probs")
	s <- 1-pp[,1]
    qtile <- unique(quantile(s, seq(0,1,by=.1)))
    if(length(qtile) == 11){
        g_fac <- cut(s, breaks=quantile(s, seq(0,1,by=.1)), right=FALSE, include.lowest=FALSE)
        w <- lapply(1:length(levels(g_fac)), function(x)which(g_fac == levels(g_fac)[x]))
        obs <- sapply(w, function(x)table(factor(as.numeric(y[x]), levels=1:length(levels(y)), labels=levels(y))))
        exp <- sapply(w, function(x)colSums(pp[x, ]))
        lt5 <- mean(exp < 5)
        if(lt5 != 0 ){warning(paste(round(lt5*100), "% of expected counts < 5", sep=""))}
        Cg <- sum(c((obs-exp)^2/exp))
        Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1)
        Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=FALSE)
    }
    else{
        Cg <- Cgrefdf <- Cg_p <- NA
    }
    predcat <- predict(obj, type="class")
	r2_count <- mean(predcat == y)
	modal <- max(table(y))
	r2_counta <- (sum(y == predcat) - modal)/(length(y) - modal)
	r2_mcf <- 1-(logLik(obj)/logLik(update(obj, .~1)))
	r2_mcfa <- 1-((logLik(obj) - (obj$edf + ncol(pp)))/logLik(update(obj, .~1)))
	g2 <- -2*(logLik(update(obj, .~1)) - logLik(obj))
	r2_ml <- 1-exp(-g2/length(y))
    r2cu <- (r2_ml)/(1-exp(2*logLik(update(obj, .~1))/nrow(pp)))
	res <- matrix(nrow=7, ncol=2)
	colnames(res) <- c("Estimate", "p-value")
	rownames(res) <- c("Fagerland, Hosmer and Bonfi", "Count R2", "Count R2 (Adj)",
		"ML R2", "McFadden R2", "McFadden R2 (Adj)", "Cragg-Uhler(Nagelkerke) R2")
	res[1,1] <- Cg
	res[1,2] <- Cg_p
	res[2,1] <- r2_count
	res[3,1] <- r2_counta
	res[4,1] <- r2_ml
	res[5,1] <- r2_mcf
	res[6,1] <- r2_mcfa
	res[7,1] <- r2cu

	permres <- NULL
	if(permute & length(qtile == 11)){
		X <- model.matrix(obj)
		l <- levels(y)
		for(i in 1:length(l)){
		y <- model.response(model.frame(obj))
		y <- relevel(y, l[i])
		tmpmod <- multinom(y ~ X-1, trace=FALSE)
		pp <- predict(tmpmod, type="probs")
		s <- 1-pp[,1]
		g_fac <- cut(s, breaks=quantile(s, seq(0,1,by=.1)), right=FALSE, include.lowest=FALSE)
		w <- lapply(1:length(levels(g_fac)), function(x)which(g_fac == levels(g_fac)[x]))
		obs <- sapply(w, function(x)table(factor(as.numeric(y[x]), levels=1:length(levels(y)), labels=levels(y))))
		exp <- sapply(w, function(x)colSums(pp[x, ]))
		Cg <- sum(c((obs-exp)^2/exp))
		Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1)
		Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=FALSE)
		permres <- rbind(permres, data.frame(base = l[i], Cg = Cg, p = Cg_p, stringsAsFactors=TRUE))
		}
	}

	if(is.null(permres)){
		out <- list(result=res)
	}
	else{
		out <- list(result=res, permres=permres)
	}
	class(out) <- "mnlfit"
	return(out)
}

##' @method print mnlfit
print.mnlfit <- function(x, ..., digits=3){
	fmt <- paste("%.", digits, "f", sep="")
	est <- sprintf(fmt, x$result[,1])
	p <- gsub("NA", "", sprintf(fmt, x$result[,2]))
	newx <- cbind(est, p)
	dimnames(newx) <- dimnames(x$result)
	cat("Fit Statistics\n")
	print(newx, quote=FALSE)
	if(exists("permres", x)){
		permbase <- as.character(x$permres[,1])
		permest <- sprintf(fmt, x$permres[,2])
		permp <- gsub("NA", "", sprintf(fmt, x$permres[,3]))
		newp <- cbind(permbase, permest, permp)
		dimnames(newp) <- dimnames(x$permres)
		cat("\nPermutations for Fagerland et. al\n")
		print(newp, quote=FALSE)
	}
}

##' @method print ordfit
print.ordfit <- function(x,..., digits=3){
	fmt <- paste("%.", digits, "f", sep="")
	est <- sprintf(fmt, x[,1])
	p <- gsub("NA", "", sprintf(fmt, x[,2]))
	newx <- cbind(est, p)
	dimnames(newx) <- dimnames(x)
	print(newx[-(1:4), 1, drop=FALSE], quote=FALSE)
}



#' Fit Statistics for Proportional Odds Logistic Regression Models
#' 
#' For objects of class \code{polr}, it calculates a number of fit statistics
#' and specification tests.
#' 
#' 
#' @aliases ordfit print.ordfit
#' @param obj A model object of class \code{polr}.
#' @return An object of class \code{ordfit} which is a matrix containing
#' statistics and specification tests.
#' @author Dave Armstrong
#' @references Lipsitz, S. R., Fitzmaurice, G. M. and Mohlenberghs, G.  1996.
#' Goodness-of-fit Tests for Ordinal Response Regression Models.  Applied
#' Statistics, 45: 175-190.\cr Pulkstenis, E. and Robinson, T. J. 2004.
#' Goodness-of-fit Test for Ordinal Response Regression Models.  Statistics in
#' Medicine, 23: 999-1014. \cr Fagerland, M. W. and Hosmer, D. W.  2013.  A
#' Goodness-of-fit Test for the Proportional Odds Regression Model.  Statistics
#' in Medicine 32(13): 2235-2249.
#' 
#' @export
#' 
#' @examples
#' 
#' library(MASS)
#' data(france)
#' polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
#' ordfit(polr.mod)
#' 
ordfit <- function(obj){
	# combfun <- function(mytab){
	# 	lt <- sum(mytab[,2] < 5)
	# 	while(lt > 0){
	# 		myw <- which(mytab[,2] < 5)
	# 		combrow <- ifelse(max(myw) == 1, 2, max(myw)-1)
	# 		mytab[combrow, ] <- apply(mytab[c(combrow, max(myw)), ], 2, sum)
	# 		mytab <- matrix(mytab[-max(myw), ], ncol=2)
	# 		lt <- sum(mytab[,2] < 5)
	# 	}
	# 	mytab
	# }
	# combcount <- function(mytab){
	# 	k <- 0
	# 	lt <- sum(mytab[,2] < 5)
	# 	while(lt > 0){
	# 		myw <- which(mytab[,2] < 5)
	# 		combrow <- ifelse(max(myw) == 1, 2, max(myw)-1)
	# 		mytab[combrow, ] <- apply(mytab[c(combrow, max(myw)), ], 2, sum)
	# 		mytab <- matrix(mytab[-max(myw), ], ncol=2)
	# 		lt <- sum(mytab[,2] < 5)
	# 		k <- k+1
	# 	}
	# 	k
	# }
	pp <- predict(obj, type="probs")
	# ints <- 1:ncol(pp)
	# n <- nrow(pp)
	# up <- floor(n/(5*ncol(pp)))
	# g <- ifelse(up > 10, 10, max(6, up))
	# s <- pp %*% ints
	# g_fac <- cut(s, breaks=(g))
	X <- model.matrix(obj)[,-1]
	y <- model.response(model.frame(obj))
	orig <- polr(y ~ X)
	# upmod <- polr(y ~ X + g_fac)
	# lrt <- lrtest(orig, upmod)
	# facs <- which(attr(obj$terms, "dataClasses") == "factor")[-1]
	# mf <- model.frame(obj)
	# pats <- factor(apply(mf[, facs, drop=F], 1, function(x)paste(x, collapse=":")))
	predcat <- predict(obj, type="class")
	# prchi2 <- 0
	# prD <- 0
	# combcountPR <- 0
	# for(i in 1:length(levels(pats))){
	# 	w <- which(pats == levels(pats)[i])
	# 	tmpg <- cut(s[w], breaks=2)
	# 	tmpdat <- data.frame(obs = y[w], expect = predcat[w], stringsAsFactors=TRUE)
	# 	m1 <- apply(tmpdat[which(tmpg == levels(tmpg)[1]), ], 2, function(x)table(factor(x, levels=levels(y))))
	# 	combcountPR <- combcountPR + combcount(m1)
	# 	lt51 <- sum(m1[,2] < 5)
	# 	while(lt51 > 0 & nrow(m1) > 1){
	# 		w1 <- which(m1[,2] < 5)
	# 		combrow <- ifelse(max(w1) == 1, 2, max(w1)-1)
	# 		m1[combrow, ] <- apply(m1[c(combrow, max(w1)), ], 2, sum)
	# 		m1 <- m1[-max(w1), ]
	# 		if(is.null(nrow(m1))){m1 <- matrix(m1, ncol=2)}
	# 		lt51 <- sum(m1[,2] < 5)
	# 	}
	# 	m2 <- apply(tmpdat[which(tmpg == levels(tmpg)[2]), ], 2, function(x)table(factor(x, levels=levels(y))))
	# 	combcountPR <- combcountPR + combcount(m2)
	# 	lt52 <- sum(m2[,2] < 5)
	# 	while(lt52 > 0 & nrow(m2) >1){
	# 		w2 <- which(m2[,2] < 5)
	# 		combrow <- ifelse(max(w2) == 1, 2, max(w2)-1)
	# 		m2[combrow, ] <- apply(m2[c(combrow, max(w2)), ], 2, sum)
	# 		m2 <- m2[-max(w2), ]
	# 		if(is.null(nrow(m2))){m2 <- matrix(m2, ncol=2)}
	# 		lt52 <- sum(m2[,2] < 5)
	# 	}
	# 	if(combcountPR > 0){warning(paste(combcountPR, " rows combined due to low expected counts for P&R tests", sep=""), call.=FALSE)}
	# 	d1 <- (m1[,1] - m1[,2])^2/m1[,2]
	# 	d2 <- (m2[,1] - m2[,2])^2/m2[,2]
	# 	prchi2 <- prchi2 + sum(d1, d2)
	# 	d1a <- m1[,1] * log(m1[,1]/m1[,2])
	# 	d2a <- m2[,1] * log(m2[,1]/m2[,2])
	# 	prD <- prD + sum(d1a, d2a)
	# }
	# prD <- 2*prD
	# refdf <- (2*nlevels(pats) -1)*(ncol(pp)-1) - length(facs) - 1
	# prchi2_p <- pchisq(prchi2, refdf, lower.tail=F)
	# prD_p <- pchisq(prD, refdf, lower.tail=F)
	# tmp <- data.frame(y = y, exp = predcat, g = g_fac, stringsAsFactors=TRUE)
	# tabs <- by(tmp[,c("y", "exp")], list(tmp$g), apply, 2, function(x)table(factor(x, levels=levels(y))))
	#
	# combk <- sum(sapply(tabs, combcount))
	# if(combk > 0){warning(paste(combk, " rows combined due to low expected counts for F&H tests", sep=""), call.=FALSE)}
	# tabs <- lapply(tabs, combfun)
	# Cg <- sum(unlist(lapply(tabs, function(x)sum((x[,1]-x[,2])^2/x[,2]))))
	# Cgrefdf <- (nlevels(g_fac)-2)*(ncol(pp)-1) + (ncol(pp) - 2)
	# Cg_p <- pchisq(Cg, Cgrefdf, lower.tail=F)
	r2_count <- mean(predcat == y)
	modal <- max(table(y))
	r2_counta <- (sum(y == predcat) - modal)/(length(y) - modal)
	b <- matrix(c(0, obj$coef), ncol=1)
	X <- model.matrix(obj)
	r2_mz <- (t(b)%*%var(X)%*%b)/(t(b)%*%var(X)%*%b + (pi^2/3))
	r2_mcf <- 1-(logLik(obj)/logLik(update(obj, .~1)))
	r2_mcfa <- 1-((logLik(obj) - obj$edf)/logLik(update(obj, .~1)))
	g2 <- -2*(logLik(update(obj, .~1)) - logLik(obj))
	r2_ml <- 1-exp(-g2/length(y))

	res <- matrix(nrow=10, ncol=2)
	colnames(res) <- c("Estimate", "p-value")
	rownames(res) <- c("Lipsitz et al", "Pulkstein and Robinson (chi-squared)", "Pulkstein and Robinson (D)", "Fagerland and Hosmer", "Count R2", "Count R2 (Adj)",
		"ML R2", "McFadden R2", "McFadden R2 (Adj)", "McKelvey & Zavoina R2")
	# res[1,1] <- lrt[2,4]
	# res[1,2] <- lrt[2,5]
	# res[2,1] <- prchi2
	# res[2,2] <- prchi2_p
	# res[3,1] <- prD
	# res[3,2] <- prD_p
	# res[4,1] <- Cg
	# res[4,2] <- Cg_p
	res[5,1] <- r2_count
	res[6,1] <- r2_counta
	res[7,1] <- r2_ml
	res[8,1] <- r2_mcf
	res[9,1] <- r2_mcfa
	res[10,1] <- r2_mz
	class(res) <- c("ordfit", "matrix")
	print(res)
}



#' Plot Average Effects of Variables in Proportional Odds Logistic Regression
#' 
#' For objects of class \code{polr} the function plots the average effect of a
#' single variable holding all other variables at their observed values.
#' 
#' Following the advice of Hanmer and Kalkan (2013) the function calculates the
#' average effect of a variable holding all other variables at observed values
#' and then plots the result.
#' 
#' @param obj An object of class \code{polr}
#' @param varname A string providing the name of the variable for which you
#' want the plot to be drawn.
#' @param data Data used to estimate \code{obj}.
#' @param R Number of simulations to generate confidence intervals.
#' @param nvals Number of evaluation points of the function
#' @param plot Logical indicating whether or not the result should be plotted
#' (if \code{TRUE}) or returned to the console (if \code{FALSE}).
#' @param returnInd Logical indicating whether average individual probabilities
#' should be returned.
#' @param returnMprob Logical indicating whether marginal probabilities,
#' averaged over individuals, should be returned.
#' @param \dots Arguments passed down to the call to \code{xyplot}
#' @return Either a plot or a list with a data frame containing the variables
#' \item{mean}{The average effect (i.e., predicted probability)}
#' \item{lower}{The lower 95\% confidence bound} \item{upper}{The upper 95\%
#' confidence bound} \item{y}{The values of the dependent variable being
#' predicted} \item{x}{The values of the independent variable being
#' manipulated} and the elements Ind or Mprob, as requested.
#' @author Dave Armstrong
#' @references Hanmer, M.J. and K.O. Kalkan.  2013.  \sQuote{Behind the Curve:
#' Clarifying the Best Approach to Calculating Predicted Probabilities and
#' Marginal Effects from Limited Dependent Variable Models}.  American Journal
#' of Political Science.  57(1): 263-277.
#' 
#' @export
#' 
#' @examples
#' 
#' library(MASS)
#' data(france)
#' polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
#' \dontrun{ordAveEffPlot(polr.mod, "lrself", data=france)	}
#' 
ordAveEffPlot <- function(obj, varname, data, R=1500, nvals=25, plot=TRUE, returnInd=FALSE, returnMprob=FALSE,...){
    pfun <- switch(obj$method, logistic = plogis, probit = pnorm,
           loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)
    rI <- rMP <- list()
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, c(coef(obj), obj$zeta), vcov(obj))
    ints <- list()
    nc <- length(obj$coef)
    for(i in 1:(length(obj$lev)-1)){
        ints[[i]] <- matrix(b[,(nc+i)], byrow=TRUE, nrow=nrow(model.matrix(obj)), ncol=R)
    }
    if(is.numeric(data[[varname]])){
	  s <- seq(min(data[[varname]], na.rm=TRUE), max(data[[varname]], na.rm=TRUE), length=nvals)
	}
    else{
        s <- obj$xlevels[[varname]]
        nvals <- length(s)
    }
    d0 <- data
    m <- l <- u <- matrix(NA, nrow=nvals, ncol=length(obj$lev))
    for(i in 1:length(s)){
        if(is.numeric(data[[varname]])){
    		d0[[varname]] <- s[i]
        }
        else{
            d0[[varname]] <- factor(i, levels=1:length(s), labels=s)
        }
       d0[[varname]][which(is.na(data[[varname]]))] <- NA
       X <- model.matrix(formula(obj), data=d0)[,-1]
       XB <- X %*% t(b[,1:length(obj$coef)])
       Q <- lapply(ints, function(x)pfun(x-XB))
       P <- list()
       for(j in 1:length(Q)){
           if(j == 1){
               P[[j]] <- Q[[j]]
           }
           else{
               P[[j]] <- Q[[j]]-Q[[(j-1)]]
           }
        }
        P[[(length(Q)+1)]] <- 1-Q[[length(Q)]]
        if(returnInd){
            rI[[i]] <- sapply(P, rowMeans)
        }
        if(returnMprob){
            rMP[[i]] <- sapply(P, colMeans)
        }
        P2 <- sapply(P, colMeans)
        m[i,] <- colMeans(P2)
        l[i,] <- apply(P2, 2, quantile, .025)
        u[i,] <- apply(P2, 2, quantile, .975)
    }

	tmp <- data.frame(
		mean = c(m),
		lower = c(l),
		upper = c(u),
		y = rep(obj$lev, each = length(s)), 
		stringsAsFactors=TRUE)
	if(is.numeric(data[[varname]])){tmp$s <- s}
	else{tmp$s <- factor(1:length(s), labels=s)}
		if(plot){
			if(is.factor(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				scales=list(x=list(at=1:length(s), labels=s)),
				panel = function(x,y, lower, upper, subscripts){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower[subscripts], x, upper[subscripts], lty=1, col="black")
			})
			}
			if(is.numeric(data[[varname]])){
			pl <- xyplot(mean ~ s | y, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				panel=panel.ci, zl=FALSE)
			}
			return(pl)
		}
		else{
            ret <- list()
            ret$data <- tmp
            if(returnInd){
                ret$Ind <- rI
            }
            if(returnMprob){
                ret$Mprob <- rMP
            }
			return(ret)
		}
	}



#' Deviance and Chi-squared Goodness-of-Fit Test for Poisson Models
#' 
#' Deviance and Chi-squared goodness-of-fit test of the null hypothesis that
#' poisson variance is appropriate to model the conditional dispersion of the
#' data, given a particular model.
#' 
#' 
#' @param obj A model object of class \code{glm} (with \code{family=poisson}).
#' @return A 2x2 data frame with rows representing the different types of
#' statistics (Deviance and Chi-squared) and columns representing the test
#' statistic and p-value.
#' @author Dave Armstrong
#' @references Dobson, A. J. (1990) An Introduction to Generalized Linear
#' Models. London: Chapman and Hall.
#' 
#' @export
#' 
#' @examples
#' 
#' ## Example taken from MASS help file for glm, identified to be
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts, stringsAsFactors=TRUE))
#' glm.D93 <- glm(counts ~ outcome + treatment, family=poisson())
#' poisGOF(glm.D93)
#' 
poisGOF <-
function (obj)
{
    if (!("glm" %in% class(obj))) {
        stop("poisGOF only works on objects of class glm (with a poisson family)\n")
    }
    if (!("y" %in% names(obj))) {
        obj <- update(obj, y = TRUE)
    }
    ind.chisq <- (((obj$y - obj$fitted)^2)/obj$fitted)
    dev <- obj$deviance
    df <- obj$df.residual
    p.chisq <- pchisq(sum(ind.chisq), df = df, lower.tail = FALSE)
    p.dev <- pchisq(dev, df = df, lower.tail = FALSE)
    vec <- sprintf("%.3f", c(sum(ind.chisq), dev, p.chisq, p.dev))
    mat <- matrix(vec, ncol = 2, nrow = 2)
    rownames(mat) <- c("Chi-squared", "Deviance")
    colnames(mat) <- c("Stat", "p-value")
    mat <- as.data.frame(mat)
    return(mat)
}


#' Proportional and Expected Proportional Reductions in Error
#' 
#' Calculates proportional reduction in error (PRE) and expected proportional
#' reduction in error (epre) from Herron (1999).
#' 
#' Proportional reduction in error is calculated as a function of correct and
#' incorrect predictions (and the probabilities of correct and incorrect
#' predictions for ePRE).  When \code{sim=TRUE}, a parametric bootstrap will be
#' used that draws from the multivariate normal distribution centered at the
#' coefficient estimates from the model and using the estimated
#' variance-covariance matrix of the estimators as Sigma.  This matrix is used
#' to form \code{R} versions of XB and predictions are made for each of the
#' \code{R} different versions of XB.  Confidence intervals can then be created
#' from the bootstrap sampled (e)PRE values.
#' 
#' @param mod1 A model of class \code{glm} (with family \code{binomial}),
#' \code{polr} or \code{multinom} for which (e)PRE will be calculated.
#' @param mod2 A model of the same class as \code{mod1} against which
#' proportional reduction in error will be measured.  If \code{NULL}, the null
#' model will be used.
#' @param sim A logical argument indicating whether a parametric bootstrap
#' should be used to calculate confidence bounds for (e)PRE.  See
#' \code{Details} for more information.
#' @param R Number of bootstrap samples to be drawn if \code{sim=TRUE}.
#' @return An object of class \code{pre}, which is a list with the following
#' elements: \item{pre}{The proportional reduction in error} \item{epre}{The
#' expected proportional reduction in error} \item{m1form}{The formula for
#' model 1} \item{m2form}{The formula for model 2} \item{pcp}{The percent
#' correctly predicted by model 1} \item{pmc}{The percent correctly predicted
#' by model 2} \item{epcp}{The expected percent correctly predicted by model 1}
#' \item{epmc}{The expected percent correctly predicted by model 2}
#' \item{pre.sim}{A vector of bootstrapped PRE values if \code{sim=TRUE}}
#' \item{epre.sim}{A vector of bootstrapped ePRE values if \code{sim=TRUE}}
#' @author Dave Armstrong
#' @references Herron, M.  1999.  Postestimation Uncertainty in Limited
#' Dependent Variable Models.  Political Analysis 8(1): 83--98.
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2), data=france, family=binomial)
#' pre(left.mod)
#' 
pre <-
function (mod1, mod2 = NULL, sim = FALSE, R = 2500)
{
    if (!is.null(mod2)) {
        if (mean(class(mod1) == class(mod2)) != 1) {
            stop("Model 2 must be either NULL or of the same class as Model 1\n")
        }
    }
    if (!any(class(mod1) %in% c("polr", "multinom", "glm"))) {
        stop("pre only works on models of class glm (with binomial family), polr or multinom\n")
    }
    if ("glm" %in% class(mod1)) {
        if (!(family(mod1)$link %in% c("logit", "probit", "cloglog",
            "cauchit"))) {
            stop("PRE only calculated for models with logit, probit, cloglog or cauchit links\n")
        }
        if (is.null(mod2)) {
            y <- mod1[["y"]]
            mod2 <- update(mod1, ". ~ 1", data=model.frame(mod1))
        }
        pred.mod2 <- as.numeric(predict(mod2, type = "response") >=
            0.5)
        pmc <- mean(mod2$y == pred.mod2)
        pred.y <- as.numeric(predict(mod1, type = "response") >=
            0.5)
        pcp <- mean(pred.y == mod1$y)
        pre <- (pcp - pmc)/(1 - pmc)
        pred.prob1 <- predict(mod1, type = "response")
        pred.prob2 <- predict(mod2, type = "response")
        epcp <- (1/length(pred.prob1)) * (sum(pred.prob1[which(mod1$y ==
            1)]) + sum(1 - pred.prob1[which(mod1$y == 0)]))
        epmc <- (1/length(pred.prob2)) * (sum(pred.prob2[which(mod2$y ==
            1)]) + sum(1 - pred.prob2[which(mod2$y == 0)]))
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1.sim <- mvrnorm(R, coef(mod1), vcov(mod1))
            b2.sim <- mvrnorm(R, coef(mod2), vcov(mod2))
            mod1.probs <- family(mod1)$linkinv(model.matrix(mod1) %*%
                t(b1.sim))
            mod2.probs <- family(mod2)$linkinv(model.matrix(mod2) %*%
                t(b2.sim))
            pmcs <- apply(mod2.probs, 2, function(x) mean(as.numeric(x >
                0.5) == mod2$y))
            pcps <- apply(mod1.probs, 2, function(x) mean(as.numeric(x >
                0.5) == mod1$y))
            pre.sim <- (pcps - pmcs)/(1 - pmcs)
            epmc.sim <- apply(mod2.probs, 2, function(x) (1/length(x)) *
                (sum(x[which(mod2$y == 1)]) + sum(1 - x[which(mod2$y ==
                  0)])))
            epcp.sim <- apply(mod1.probs, 2, function(x) (1/length(x)) *
                (sum(x[which(mod1$y == 1)]) + sum(1 - x[which(mod1$y ==
                  0)])))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    if ("multinom" %in% class(mod1)) {
        mod1 <- update(mod1, ".~.", model = TRUE, trace = FALSE)
        if (is.null(mod2)) {
            mod2 <- update(mod1, ". ~ 1", data = mod1$model,
                trace = FALSE)
        }
        pred.prob1 <- predict(mod1, type = "prob")
        pred.prob2 <- predict(mod2, type = "prob")
        pred.cat1 <- apply(pred.prob1, 1, which.max)
        pred.cat2 <- apply(pred.prob2, 1, which.max)
        pcp <- mean(as.numeric(mod1$model[, 1]) == pred.cat1)
        pmc <- max(table(as.numeric(mod2$model[, 1]))/sum(table(mod2$model[,
            1])))
        pre <- (pcp - pmc)/(1 - pmc)
        epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), as.numeric(mod1$model[,
            1]))])
        tab <- table(mod1$model[, 1])/sum(table(mod1$model[,
            1]))
        epmc <- mean(tab[as.numeric(mod2$model[, 1])])
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1.sim <- mvrnorm(R, c(t(coef(mod1))), vcov(mod1))
            b2.sim <- mvrnorm(R, c(t(coef(mod2))), vcov(mod2))
            mod.levs <- rownames(coef(mod1))
            var.levs <- rownames(contrasts(mod1$model[, 1]))
            tmp.reord <- c(match(mod.levs, var.levs), (1:length(var.levs))[-match(mod.levs,
                var.levs)])
            reord <- match(1:length(var.levs), tmp.reord)
            tmp <- lapply(1:nrow(b1.sim), function(x) matrix(b1.sim[x,
                ], ncol = ncol(coef(mod1)), byrow = TRUE))
            tmp2 <- lapply(tmp, function(x) cbind(model.matrix(mod1) %*%
                t(x), 0))
            tmp2 <- lapply(tmp2, function(x) x[, reord])
            mod1.probs <- lapply(tmp2, function(x) exp(x)/apply(exp(x),
                1, sum))
            tmp <- lapply(1:nrow(b2.sim), function(x) matrix(b2.sim[x,
                ], ncol = ncol(coef(mod2)), byrow = TRUE))
            tmp2 <- lapply(tmp, function(x) cbind(model.matrix(mod2) %*%
                t(x), 0))
            tmp2 <- lapply(tmp2, function(x) x[, reord])
            mod2.probs <- lapply(tmp2, function(x) exp(x)/apply(exp(x),
                1, sum))
            pred.cat1 <- lapply(mod1.probs, function(x) apply(x,
                1, which.max))
            pred.cat2 <- lapply(mod2.probs, function(x) apply(x,
                1, which.max))
            pcp.sim <- sapply(pred.cat1, function(x) mean(x ==
                as.numeric(mod1$model[, 1])))
            pmc.sim <- sapply(pred.cat2, function(x) mean(x ==
                as.numeric(mod1$model[, 1])))
            pre.sim <- (pcp.sim - pmc.sim)/(1 - pmc.sim)
            epcp.sim <- sapply(mod1.probs, function(z) mean(z[cbind(1:nrow(mod1$model),
                as.numeric(mod1$model[, 1]))]))
            epmc.sim <- sapply(mod2.probs, function(z) mean(z[cbind(1:nrow(mod1$model),
                as.numeric(mod1$model[, 1]))]))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    if ("polr" %in% class(mod1)) {
        if (is.null(mod1$Hess)) {
            mod1 <- update(mod1, Hess = TRUE)
        }
        if (is.null(mod2)) {
            mod2 <- update(mod1, ". ~ 1", data = mod1$model,
                model = TRUE, Hess = TRUE)
        }
        pred.prob1 <- predict(mod1, type = "prob")
        pred.prob2 <- predict(mod2, type = "prob")
        pred.cat1 <- apply(pred.prob1, 1, which.max)
        pred.cat2 <- apply(pred.prob2, 1, which.max)
        pcp <- mean(as.numeric(mod1$model[, 1]) == pred.cat1)
        pmc <- max(table(as.numeric(mod2$model[, 1]))/sum(table(mod2$model[,
            1])))
        pre <- (pcp - pmc)/(1 - pmc)
        epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), as.numeric(mod1$model[,
            1]))])
        tab <- table(mod1$model[, 1])/sum(table(mod1$model[,
            1]))
        epmc <- mean(tab[as.numeric(mod2$model[, 1])])
        epre <- (epcp - epmc)/(1 - epmc)
        if (sim) {
            b1 <- c(coef(mod1), mod1$zeta)
            b2 <- c(coef(mod2), mod2$zeta)
            v1 <- invisible(vcov(mod1))
            v2 <- invisible(vcov(mod2))
            b1.sim <- mvrnorm(R, b1, v1)
            b2.sim <- mvrnorm(R, b2, v2)
            mod1.probs <- lapply(1:nrow(b1.sim), function(x) simPredpolr(mod1,
                b1.sim[x, ], n.coef = length(coef(mod1))))
            mod2.probs <- lapply(1:nrow(b2.sim), function(x) simPredpolr(mod2,
                b2.sim[x, ], n.coef = length(coef(mod2))))
            pred.cat1 <- lapply(mod1.probs, function(x) apply(x,
                1, which.max))
            pred.cat2 <- lapply(mod2.probs, function(x) apply(x,
                1, which.max))
            pcp.sim <- sapply(pred.cat1, function(x) mean(x ==
                as.numeric(mod1$model[, 1])))
            pmc.sim <- sapply(pred.cat2, function(x) mean(x ==
                as.numeric(mod1$model[, 1])))
            pre.sim <- (pcp.sim - pmc.sim)/(1 - pmc.sim)
            epcp.sim <- sapply(mod1.probs, function(z) mean(z[cbind(1:nrow(mod1$model),
                as.numeric(mod1$model[, 1]))]))
            epmc.sim <- sapply(mod2.probs, function(z) mean(z[cbind(1:nrow(mod1$model),
                as.numeric(mod1$model[, 1]))]))
            epre.sim <- (epcp.sim - epmc.sim)/(1 - epmc.sim)
        }
    }
    ret <- list()
    ret$pre <- pre
    ret$epre <- epre
    form1 <- formula(mod1)
    form2 <- formula(mod2)
    ret$m1form <- paste(form1[2], form1[1], form1[3], sep = " ")
    ret$m2form <- paste(form2[2], form2[1], form2[3], sep = " ")
    ret$pcp <- pcp
    ret$pmc <- pmc
    ret$epmc <- epmc
    ret$epcp <- epcp
    if (sim) {
        ret$pre.sim <- pre.sim
        ret$epre.sim <- epre.sim
    }
    class(ret) <- "pre"
    return(ret)
}


#' Print method for objects of class pre
#' 
#' Prints the output from an object of class \code{pre}.  The function prints
#' all components of the calculation and optionally simulated confidence
#' bounds.
#' 
#' 
#' @param x An object of class \code{pre}.
#' @param sim.ci Coverage for the simulated confidence interval, if
#' \code{sim=TRUE} in the call to \code{pre}.
#' @param ... Other arguments passed to print, currently not implemented
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @method print pre
#' @seealso \code{pre}
print.pre <-
function (x, ..., sim.ci = 0.95)
{
    cat("mod1: ", as.character(x$m1form), "\n")
    cat("mod2: ", as.character(x$m2form), "\n\n")
    cat("Analytical Results\n")
    cat(" PMC = ", sprintf("%2.3f", x$pmc), "\n")
    cat(" PCP = ", sprintf("%2.3f", x$pcp), "\n")
    cat(" PRE = ", sprintf("%2.3f", x$pre), "\n")
    cat("ePMC = ", sprintf("%2.3f", x$epmc), "\n")
    cat("ePCP = ", sprintf("%2.3f", x$epcp), "\n")
    cat("ePRE = ", sprintf("%2.3f", x$epre), "\n\n")
    low <- (1 - sim.ci)/2
    up <- 1 - low
    if (exists("pre.sim", x)) {
        pre.ci <- sprintf("%2.3f", quantile(x$pre.sim, c(0.5,
            low, up)))
        epre.ci <- sprintf("%2.3f", quantile(x$epre.sim, c(0.5,
            low, up)))
        tmp <- rbind(pre.ci, epre.ci)
        rownames(tmp) <- c(" PRE", "ePRE")
        colnames(tmp) <- c("median", "lower", "upper")
        cat("Simulated Results\n")
        print(tmp, quote = FALSE)
    }
}
#' 
#' @export
#' 
probit_cc <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- X %*% b
    phat <- pnorm(xb)
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    probitcc <- (b[int.var] - (b[vars[1]] + b[int.var] * X[,
        vars[2]]) * (b[vars[2]] + b[int.var] * X[, vars[1]]) *
        xb) * phi
    d2f <- -xb * phi
    d3f <- (xb^2 - 1) * phi
    b1b4x2 <- b[vars[1]] + b[int.var] * X[, vars[2]]
    b2b4x1 <- b[vars[2]] + b[int.var] * X[, vars[1]]
    deriv11 <- b[int.var] * d2f * X[, vars[1]] + b2b4x1 * d2f +
        b1b4x2 * b2b4x1 * d3f * X[, vars[1]]
    deriv22 <- b[int.var] * d2f * X[, vars[2]] + b1b4x2 * d2f +
        b1b4x2 * b2b4x1 * X[, vars[2]] * d3f
    deriv44 <- phi + b[int.var] * d2f * X[, vars[1]] * X[, vars[2]] +
        X[, vars[2]] * b2b4x1 * d2f + X[, vars[1]] * b1b4x2 *
        d2f + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] *
        d3f
    derivcc <- b[int.var] * d2f + b1b4x2 * b2b4x1 * d3f
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) b[int.var] * d2f * x +
        b1b4x2 * b2b4x1 * x * d3f)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv11, deriv22, deriv44, nn, derivcc)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitcc/probit_se
    out <- data.frame(int_eff = probitcc, linear = linear, phat = phat,
        se_int_eff = probit_se, zstat = probit_t, stringsAsFactors=TRUE)
    invisible(out)
}
#' 
#' @export
#' 
probit_cd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- X %*% b
    phat <- pnorm(xb)
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    dum <- vars[which(sapply(apply(X[, vars], 2, table), length) ==
        2)]
    cont <- vars[which(vars != dum)]
    X1 <- X2 <- X
    X1[, dum] <- 1
    X1[, int.var] <- X1[, cont] * X1[, dum]
    phi1 <- dnorm(X1 %*% b)
    d2f1 <- -(X1 %*% b) * phi1
    ie1 <- (b[cont] + b[int.var]) * phi1
    X2[, dum] <- 0
    X2[, int.var] <- X2[, cont] * X2[, dum]
    phi2 <- dnorm(X2 %*% b)
    d2f2 <- -(X2 %*% b) * phi2
    ie2 <- b[cont] * phi2
    probitcd <- ie1 - ie2
    deriv1 <- phi1 - phi2 + b[cont] * X[, cont] * (d2f1 - d2f2) +
        b[int.var] * X[, cont] * d2f1
    deriv2 <- (b[cont] + b[int.var]) * d2f1
    deriv3 <- phi1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
    deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) ((b[cont] + b[int.var]) *
        d2f1 - b[cont] * d2f2) * x)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitcd/probit_se
    out <- data.frame(int_eff = probitcd, linear = linear, phat = phat,
        se_int_eff = probit_se, zstat = probit_t, stringsAsFactors=TRUE)
    invisible(out)
}
#' 
#' @export
#' 
probit_dd <-
function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
{
    xb <- X %*% b
    phat <- pnorm(xb)
    phi <- dnorm(xb)
    linear <- b[int.var] * phi
    X11 <- X01 <- X10 <- X00 <- X
    X11[, vars[1]] <- 1
    X11[, vars[2]] <- 1
    X10[, vars[1]] <- 1
    X10[, vars[2]] <- 0
    X01[, vars[1]] <- 0
    X01[, vars[2]] <- 1
    X00[, vars[1]] <- 0
    X00[, vars[2]] <- 0
    X00[, int.var] <- X00[, vars[1]] * X00[, vars[2]]
    X11[, int.var] <- X11[, vars[1]] * X11[, vars[2]]
    X01[, int.var] <- X01[, vars[1]] * X01[, vars[2]]
    X10[, int.var] <- X10[, vars[1]] * X10[, vars[2]]
    Xb11 <- X11 %*% b
    Xb00 <- X00 %*% b
    Xb10 <- X10 %*% b
    Xb01 <- X01 %*% b
    phat11 <- pnorm(Xb11)
    phat00 <- pnorm(Xb00)
    phat10 <- pnorm(Xb10)
    phat01 <- pnorm(Xb01)
    phi11 <- dnorm(Xb11)
    phi00 <- dnorm(Xb00)
    phi10 <- dnorm(Xb10)
    phi01 <- dnorm(Xb01)
    probitdd <- (phat11 - phat10) - (phat01 - phat00)
    deriv1 <- phi11 - phi10
    deriv2 <- phi11 - phi01
    deriv3 <- phi11
    deriv0 <- (phi11 - phi01) - (phi10 - phi00)
    others <- X[, -c(1, match(c(vars, int.var), names(b)))]
    if (!("matrix" %in% class(others))) {
        others <- matrix(others, nrow = nrow(X))
    }
    colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
        names(b)))]
    nn <- apply(others, 2, function(x) ((phi11 - phi01) - (phi10 -
        phi00)) * x)
    nn <- array(nn, dim=dim(others))
    dimnames(nn) <- dimnames(others)
    mat123 <- cbind(deriv1, deriv2, deriv3, nn, deriv0)[,,drop=F]
    colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
    mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F]
    probit_se <- sqrt(diag(mat123 %*% vcov(obj) %*% t(mat123)))
    probit_t <- probitdd/probit_se
    out <- data.frame(int_eff = probitdd, linear = linear, phat = phat,
        se_int_eff = probit_se, zstat = probit_t, stringsAsFactors=TRUE)
    invisible(out)
}


#' Search Variable Labels Attribute
#' 
#' Data imported from SPSS or Stata comes with the variable labels set (if they
#' were set in the original dataset) as one of the dataframe's attributes.
#' This allows you to search the variable labels and returns the variable
#' column number, name and label for all variables that have partially match
#' the search term either in their labels or names.
#' 
#' For an imported Stata dataset, variable labels are in the \code{var.labels}
#' attribute of the dataset and in an SPSS dataset, they are in the
#' \code{variable.labels} attribute.  These are searched, ignoring case, for
#' the desired string
#' 
#' @aliases searchVarLabels searchVarLabels.data.frame searchVarLabels.tbl_df
#' @param dat a data frame whose variable labels you want to search.
#' @param str string used to search variable labels.
#' @return \item{matrix}{A matrix of dimensions n-matches x 2 is returned,
#' where the first column is the column number of the matching variable and the
#' second column is the variable label.  The row names of the matrix are the
#' variable names.}
#' 
#' @export
#' 
#' @author Dave Armstrong
searchVarLabels <- function(dat, str) UseMethod("searchVarLabels")

#' @export
#' @method searchVarLabels data.frame
searchVarLabels.data.frame <-
function (dat, str)
{
    if ("var.labels" %in% names(attributes(dat))) {
        vlat <- "var.labels"
    }
    if ("variable.labels" %in% names(attributes(dat))) {
        vlat <- "variable.labels"
    }
    ind <- sort(union(grep(str, attr(dat, vlat), ignore.case = TRUE), grep(str, names(dat), ignore.case = TRUE)))
    vldf <- data.frame(ind = ind, label = attr(dat, vlat)[ind], stringsAsFactors=TRUE)
    rownames(vldf) <- names(dat)[ind]
    vldf
}

#' @export
#' @method searchVarLabels tbl_df
searchVarLabels.tbl_df <-
function (dat, str)
{
    vlat <- unlist(sapply(1:ncol(dat), function(i)attr(dat[[i]], "label")))
    ind <- sort(union(grep(str, vlat, ignore.case = TRUE), grep(str, names(dat), ignore.case = TRUE)))
    vldf <- data.frame(ind = ind, label = vlat[ind], stringsAsFactors=TRUE)
    rownames(vldf) <- names(dat)[ind]
    vldf
}


#' Calculate Predictions for Proportional Odds Logistic Regression
#' 
#' Calculates predicted probabilities from models of class \code{polr} from a
#' model object and a vector of coefficient values.  This is an auxiliary
#' function used in \code{pre} if \code{sim=TRUE}.
#' 
#' 
#' @param object An object of class \code{polr}.
#' @param n.coef Number of coefficients (minus intercepts) for the \code{polr}
#' model.
#' @param coefs A vector of coefficients where elements 1 to \code{n.coef} give
#' model coefficients and elements \code{n.coef}+1 to k have intercepts.
#' @return An n x m-category matrix of predicted probabilities
#' 
#' @export
#' 
#' @author Dave Armstrong
simPredpolr <-
function (object, coefs, n.coef)
{
    X <- model.matrix(object)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L)
        X <- X[, -xint, drop = FALSE]
    n <- nrow(X)
    q <- length(coefs) - n.coef
    eta <- {
        if (n.coef > 0) {
            drop(X %*% coefs[1:n.coef])
        }
        else {
            rep(0, n)
        }
    }
    pfun <- switch(object$method, logistic = plogis, probit = pnorm,
        cauchit = pcauchy)
    cumpr <- matrix(pfun(matrix(coefs[(n.coef + 1):length(coefs)],
        n, q, byrow = TRUE) - eta), q)
    Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    return(Y)
}


#' Maximal First Differences for Zero-Inflated Models
#' 
#' Calculates the change in predicted counts or optionally the predicted
#' probability of being in the zero-count group, for maximal discrete changes
#' in all covariates holding all other variables constant at typical values.
#' 
#' The function calculates the changes in predicted counts, or optionally the
#' predicted probability of being in the zero group, for maximal discrete
#' changes in the covariates.  This function works with polynomials specified
#' with the \code{poly} function.  It also works with multiplicative
#' interactions of the covariates by virtue of the fact that it holds all other
#' variables at typical values.  By default, typical values are the median for
#' quantitative variables and the mode for factors. The way the function works
#' with factors is a bit different.  The function identifies the two most
#' different levels of the factor and calculates the change in predictions for
#' a change from the level with the smallest prediction to the level with the
#' largest prediction.
#' 
#' @param obj A model object of class \code{zeroinfl}.
#' @param data Data frame used to fit \code{object}.
#' @param typical.dat Data frame with a single row containing values at which
#' to hold variables constant when calculating first differences.  These values
#' will be passed to \code{predict}, so factors must take on a single value,
#' but have all possible levels as their levels attribute.
#' @param type Character string of either \sQuote{count} (to obtain changes in
#' predicted counts) or \sQuote{zero} (to obtain changes in the predicted
#' probability of membership in the zero group).
#' @return A list with the following elements: \item{diffs}{A matrix of
#' calculated first differences} \item{minmax}{A matrix of values that were
#' used to calculate the predicted changes}
#' @author Dave Armstrong
#' 
#' @export
#' 
ziChange <-
function (obj, data, typical.dat = NULL, type = "count")
{
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    vars.type <- as.character(formula(obj))[3]
    vars.type <- gsub("*", "+", vars.type, fixed = TRUE)
    vars.type <- strsplit(vars.type, split = "|", fixed = TRUE)[[1]][ifelse(type ==
        "count", 1, 2)]
    vars.type <- c(unlist(strsplit(vars.type, "+", fixed = TRUE)))
    vars.type <- unique(trim(vars.type))
    pols <- grep("poly", vars.type)
    if (length(pols) > 0) {
        poly.split <- strsplit(vars.type[pols], split = "")
        start <- lapply(poly.split, function(x) grep("(", x,
            fixed = TRUE) + 1)
        stop <- lapply(poly.split, function(x) grep(",", x, fixed = TRUE) -
            1)
        pol.vars.type <- sapply(1:length(poly.split), function(x) paste(poly.split[[x]][start[[x]]:stop[[x]]],
            collapse = ""))
        vars.type[pols] <- pol.vars.type
    }
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    minmax <- lapply(vars, function(x) c(NA, NA))
    meds <- lapply(vars, function(x) NA)
    names(minmax) <- names(meds) <- vars
    levs <- obj$levels
    if (length(levs) > 0) {
        for (i in 1:length(levs)) {
            tmp.levs <- paste(names(levs)[i], unlist(levs[i]),
                sep = "")
            tmp.coefs <- obj$coef[[type]][match(tmp.levs, names(obj$coef[[type]]))]
            tmp.coefs[which(is.na(tmp.coefs))] <- 0
            mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
            minmax[[names(levs)[i]]] <- factor(levs[[i]][mm],
                levels = levs[[i]])
            tmp.tab <- table(data[[names(levs)[i]]])
            meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)],
                levels = levs[[i]])
        }
    }
    vars <- vars[sapply(minmax, function(x) is.na(x[1]))]
    for (i in 1:length(vars)) {
        minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm = TRUE)
        meds[[vars[i]]] <- median(data[[vars[i]]], na.rm = TRUE)
    }
    tmp.df <- do.call(data.frame, c(lapply(meds, function(x) rep(x,
        length(meds) * 2)), stringsAsFactors=TRUE))
    if (!is.null(typical.dat)) {
        notin <- which(!(names(typical.dat) %in% names(tmp.df)))
        if (length(notin) > 0) {
            cat("The following variables in typical.dat were not found in the prediction data: ",
                names(typical.dat)[notin], "\n\n", sep = "")
            typical.dat <- typical.dat[, -notin]
        }
        for (j in 1:ncol(typical.dat)) {
            tmp.df[[names(typical.dat)[j]]] <- typical.dat[1,
                j]
            meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1,
                j])
        }
    }
    inds <- seq(1, nrow(tmp.df), by = 2)
    for (j in 1:length(minmax)) {
        tmp.df[inds[j]:(inds[j] + 1), j] <- minmax[[j]]
    }
    preds <- matrix(predict(obj, newdata = tmp.df, type = type),
        ncol = 2, byrow = TRUE)
    diffs <- cbind(preds, apply(preds, 1, diff))
    colnames(diffs) <- c("min", "max", "diff")
    rownames(diffs) <- rn
    diffs <- diffs[which(rownames(diffs) %in% vars.type), ]
    minmax.mat <- do.call(cbind, minmax)
    minmax.mat <- rbind(c(unlist(meds)), minmax.mat)
    rownames(minmax.mat) <- c("typical", "min", "max")
    ret <- list(diffs = diffs, minmax = minmax.mat)
    class(ret) <- "change"
    return(ret)
}


#' Fitted Values and CIs for 2-Categorical Interactions
#' 
#' This function makes a table of fitted values and confidence intervals for
#' all of the combinations of two categorical variables in an interaction.
#' 
#' 
#' @param eff.obj An object generated by \code{effect} from the \code{effects}
#' package where the effect is calculated for two factors involved in an
#' interaction.
#' @param digits Number of digits of the fitted values and confidence intervals
#' to print.
#' @param rownames An optional vector of row names for the table, if
#' \code{NULL}, the levels of the factor will be used
#' @param colnames An optional vector of column names for the table, if
#' \code{NULL}, the levels of the factor will be used
#' @return A matrix of fitted values and confidence intervals
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' library(effects)
#' data(Duncan, package="carData")
#' Duncan$inc.cat <- cut(Duncan$income, 3)
#' mod <- lm(prestige ~ inc.cat*type + income, data=Duncan)
#' e1 <- effect("inc.cat*type", mod)
#' cat2Table(e1)
#' 
cat2Table <- function(eff.obj, digits=2, rownames = NULL, colnames=NULL){
	print.digs <- paste("%.", digits, "f", sep="")
	cis <- paste("(",
		sprintf(print.digs, eff.obj$lower),
		",",
		sprintf(print.digs, eff.obj$upper),
		")", sep="")
	cis.mat <- matrix(cis,
		nrow = length(levels(eff.obj$x[[1]])),
		ncol = length(levels(eff.obj$x[[2]])))

	out.mat <- NULL
	est.mat <- matrix(sprintf(print.digs, eff.obj$fit),
		nrow = length(levels(eff.obj$x[[1]])),
		ncol = length(levels(eff.obj$x[[2]])))
	for(i in 1:nrow(est.mat)){
		out.mat <- rbind(out.mat, est.mat[i,], cis.mat[i,])
	}
	if(is.null(colnames)){
		colnames(out.mat) <- levels(eff.obj$x[[2]])
	}
	if(is.null(rownames)){
		rn <- levels(eff.obj$x[[1]])
	}else{
		rn <- rownames
		}
		rn2 <- NULL
		for(i in 1:length(rn)){
			rn2 <- c(rn2, rn[i], paste(rep(" ", i), collapse=""))
		}
		rownames(out.mat) <- rn2
	out.mat
}


#' Tests the five Berry, Golder and Milton (2012) Interactive Hypothesis
#' 
#' This function tests the five hypotheses that Berry, Golder and Milton
#' identify as important when two quantitative variables are interacted in a
#' linear model.
#' 
#' 
#' @param obj An object of class \code{lm}.
#' @param vars A vector of two variable names giving the two quantitative
#' variables involved in the interaction.  These variables must be involved in
#' one, and only one, interaction.
#' @param digits Number of digits to be printed in the summary.
#' @param level Type I error rate for the tests.
#' @param two.sided Logical indicating whether the tests should be two-sided
#' (if \code{TRUE}, the default) or one-sided (if \code{FALSE}).
#' @return A matrix giving five t-tests.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Duncan, package="carData")
#' mod <- lm(prestige ~ income*education + type, data=Duncan)
#' BGMtest(mod, c("income", "education"))
#' 
BGMtest <- function(obj, vars, digits = 3, level = 0.05, two.sided=TRUE){
	cl <- attr(terms(obj), "dataClasses")
    if(any(cl[vars] != "numeric"))stop("Both variables in vars must be numeric")
    facs <- apply(attr(terms(obj), "factors"), 1, sum)
	if(any(facs[vars] != 2))stop("Each variable in vars must be involved in only 2 terms")
    r.x <- range(obj$model[, vars[1]])
	r.z <- range(obj$model[, vars[2]])
	b <- coef(obj)
	V <- vcov(obj)
	nb <- names(b)
	inds <- lapply(vars, function(x)grep(x, names(b)))
	A <- matrix(0, nrow=5, ncol=length(b))
    A[1:2, inds[[1]]] <- cbind(1, r.z)
	A[3:4, inds[[2]]] <- cbind(1, r.x)
    A[5, do.call(intersect, inds)] <- 1
	cond.b <- A %*% b
	sp1 <- do.call(max, lapply(strsplit(gsub("-", "", as.character(cond.b)), split=".", fixed=TRUE), function(x)nchar(x[1])))
	cond.se <- sqrt(diag(A %*% V %*% t(A)))
	sp2 <- do.call(max, lapply(strsplit(as.character(cond.se), split=".", fixed=TRUE), function(x)nchar(x[1])))
	cond.t <- cond.b/cond.se
	sp3 <- do.call(max, lapply(strsplit(gsub("-", "", as.character(cond.t)), split=".", fixed=TRUE), function(x)nchar(x[1])))
	cond.p <- (2^two.sided)*pt(abs(cond.t), obj$df.residual, lower.tail=FALSE)
	pdigs1 <- paste("%", (sp1), ".", (digits), "f", sep="")
	pdigs2 <- paste("%", (sp2), ".", (digits), "f", sep="")
	pdigs3 <- paste("%", (sp3), ".", (digits), "f", sep="")
	pdigs4 <- paste("%.", digits, "f", sep="")
	rn <- c("P(X|Zmin)", "P(X|Zmax)", "P(Z|Xmin)", "P(Z|Xmax)", "P(XZ)")
	out <- cbind(sprintf(pdigs1, cond.b), sprintf(pdigs2, cond.se),
		sprintf(pdigs3, cond.t), sprintf(pdigs4, cond.p))
	if(any(out[,1] < 0)){
		indp <- which(out[,1] > 0)
		out[indp,1] <- paste(" ", out[indp,1], sep="")
		out[indp,3] <- paste(" ", out[indp,3], sep="")
	}
	pad <- apply(out, 2, function(x)max(nchar(x)))
	nchars <- nchar(out)
	newchars <- t(apply(nchars, 1, function(x)pad-x))
	tmp <- apply(newchars, c(1,2), function(x)ifelse(x > 0, rep(" ", x), ""))
	out <- matrix(paste(tmp, out, sep=""), nrow=nrow(out), ncol=ncol(out))
	rownames(out) <- rn
	colnames(out) <- c(
		paste(paste(rep(" ", pad[1]-3), collapse=""), "est", sep=""),
		paste(paste(rep(" ", pad[2]-2), collapse=""), "se", sep=""),
		paste(paste(rep(" ", pad[3]-1), collapse=""), "t", sep=""),
		ifelse(pad[4] > 6,
			paste(paste(rep(" ", pad[4]-6), collapse=""), "p-value", collapse=""),
			"p-value"))
	print(out, quote=FALSE)
}



#' Predictions for Factor-Numeric Interactions in Linear Models
#' 
#' This function works on linear models with a single interaction between a
#' continuous (numeric) variable and a factor.  The output is a data frame that
#' gives the predicted effect of moving from each category to each other
#' category of the factor over the range of values of the continuous
#' conditioning variable.
#' 
#' 
#' @aliases intQualQuant 
#' @param obj An object of class \code{lm}.
#' @param vars A vector of two variable names giving the two quantitative
#' variables involved in the interaction.  These variables must be involved in
#' one, and only one, interaction.
#' @param level Confidence level desired for lower and upper bounds of
#' confidence interval.
#' @param varcov A potentially clustered or robust variance-covariance matrix
#' of parameters used to calculate standard errors.  If \code{NULL}, the
#' \code{vcov} function will be used.
#' @param labs An optional vector of labels that will be used to identify the
#' effects, if \code{NULL}, the factor levels will be used.
#' @param n Number of values of the conditioning variable to use.
#' @param onlySig Logical indicating whether only contrasts with significant
#' differences should be returned.  Significance is determined to exist if the
#' largest lower bound is greater than zero or the smallest upper bound is
#' smaller than zero.
#' @param type String indicating whether the conditional partial effect of the
#' factors is plotted (if \sQuote{facs}), or the conditional partial effect of
#' the quantitative variable (if \sQuote{slopes}) is produced.
#' @param plot Logical indicating whether graphical results (if \code{TRUE}) or
#' numerical results (if \code{FALSE}) are produced.
#' @param vals A vector of values at which the continuous variable will be held
#' constant.  If \code{NULL}, a sequence of length \code{n} across the
#' variable's range will be used.
#' @param rug Logical indicating whether rug plots should be plotted in the
#' panels.
#' @param ci Logical indicating whether confidence bounds should be drawn.
#' @param digits Number indicating how many decimal places to round the numeric
#' output.
#' @param ... Other arguments to be passed down to \code{effect} if
#' \code{plot.type} = \sQuote{slopes}.
#' @return For \code{type} = \sQuote{facs} and \code{plot} = \code{FALSE}, a
#' data frame with the following values: \item{fit}{The expected difference
#' between the two factor levels at the specified value of the conditioning
#' variable.} \item{se.fit}{The standard error of the expected differences. }
#' \item{x}{The value of the continuous conditioning variable}
#' \item{contrast}{A factor giving the two values of the factor being
#' evaluated.} \item{lower}{The lower 95\% confidence interval for \code{fit}}
#' \item{upper}{The upper 95\% confidence interval for \code{fit}} For
#' \code{type} = \sQuote{facs} and \code{plot} = \code{TRUE}, a lattice display
#' is returned For \code{type} = \sQuote{slopes} and \code{plot} =
#' \code{FALSE}, A character matrix with the following columns: \item{B}{The
#' conditional effect of the quantitative variable for each level of the
#' factor.} \item{SE(B)}{The standard error of the conditional effect.}
#' \item{t-stat}{The t-statistic of the conditional effect.}
#' \item{Pr(>|t|)}{The two-sided p-value.} For \code{type} = \sQuote{slopes}
#' and \code{plot} = \code{TRUE}, a lattice display is returned
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' Prestige$income <- Prestige$income/1000
#' mod <- lm(prestige ~ income * type + education, data=Prestige)
#' intQualQuant(mod, c("income", "type"), n=10, 
#' 	plot.type="none")
#' intQualQuant(mod, c("income", "type"), n=10, 
#' 	plot.type="facs")
#' intQualQuant(mod, c("income", "type"), n=10, 
#' 	plot.type="slopes")
#' 
intQualQuant <- function(obj, vars, level = .95 , varcov=NULL,
	labs = NULL, n = 10 , onlySig = FALSE, type = c("facs", "slopes"),
	plot=TRUE, vals = NULL, rug=TRUE, ci=TRUE, digits=3,...){
type=match.arg(type)
cl <- attr(terms(obj), "dataClasses")[vars]
if(length(cl) != 2){
	stop("vars must identify 2 and only 2 model terms")
}
if(!all(c("numeric", "factor") %in% cl)){
	stop("vars must have one numeric and one factor")
}
facvar <- names(cl)[which(cl == "factor")]
quantvar <- names(cl)[which(cl == "numeric")]
faclevs <- obj$xlevels[[facvar]]
if(is.null(labs)){
	labs <- faclevs
}
if(!is.null(vals)){n <- length(vals)}
if(!is.null(vals)){
	quantseq <- vals
}
else{
	qrange <- range(obj$model[[quantvar]], na.rm=TRUE)
	quantseq <- seq(qrange[1], qrange[2], length=n)
}
b <- coef(obj)
if(is.null(varcov)){
    varcov <- vcov(obj)
}
faccoef <- paste(facvar, faclevs, sep="")
main.ind <- sapply(faccoef, function(x)
	grep(paste("^", x, "$", sep=""), names(b)))
main.ind <- sapply(main.ind, function(x)
	ifelse(length(x) == 0, 0, x))
int.ind1 <- sapply(faccoef, function(x){
	g1 <- grep(paste("[a-zA-Z0-9]*\\:", x, "$", sep=""),
		names(b))
	ifelse(length(g1) == 0, 0, g1)
})
int.ind2 <- sapply(faccoef, function(x){
	g2 <- grep(paste("^", x, "\\:[a-zA-Z0-9]*", sep=""),
		names(b))
	ifelse(length(g2) == 0, 0, g2)
})
{if(sum(int.ind1) != 0){int.ind <- int.ind1}
else{int.ind <- int.ind2}}

inds <- cbind(main.ind, int.ind)
nc <- ncol(inds)
dn <- dimnames(inds)
if(!("matrix" %in% class(inds))){inds <- matrix(inds, nrow=1)}
outind <- which(main.ind == 0)
inds <- inds[-outind, ]
if(length(inds) == nc){
	inds <- matrix(inds, ncol=nc)
}
rownames(inds) <- dn[[1]][-outind]
colnames(inds) <- dn[[2]]

if(length(faclevs) < 2){stop("Factor must have at least two unique values")}
{if(length(faclevs) > 2){
combs <- combn(length(faclevs)-1, 2)
}
else{
	combs <- matrix(1:length(faclevs), ncol=1)
}}
mf <- model.frame(obj)
c2 <- combn(1:length(faclevs), 2)
dc <- dim(c2)
fl2 <- matrix(faclevs[c2], nrow=dc[1], ncol=dc[2])
l <- list()
for(i in 1:ncol(c2)){
	l[[i]] <- list()
	l[[i]][[fl2[[1,i]]]] <- mf[which(mf[[facvar]] == fl2[1,i]), quantvar]
	l[[i]][[fl2[[2,i]]]] <- mf[which(mf[[facvar]] == fl2[2,i]), quantvar]
}

tmp.A <- matrix(0, nrow=length(quantseq), ncol=length(b))
A.list <- list()
k <- 1
for(i in 1:nrow(inds)){
	A.list[[k]] <- tmp.A
	A.list[[k]][,inds[i,1]] <- 1
	A.list[[k]][,inds[i,2]] <- quantseq
	k <- k+1
}
if(nrow(inds) > 1){
for(i in 1:ncol(combs)){
	A.list[[k]] <- tmp.A
	A.list[[k]][,inds[combs[1,i], 1]] <- -1
	A.list[[k]][,inds[combs[2,i], 1]] <- 1
	A.list[[k]][,inds[combs[1,i], 2]] <- -quantseq
	A.list[[k]][,inds[combs[2,i], 2]] <- quantseq
	k <- k+1
}
}

effs <- lapply(A.list, function(x)x%*%b)
se.effs <- lapply(A.list, function(x)sqrt(diag(x %*% varcov %*%t(x))))
allcombs <- combn(length(faclevs), 2)
list.labs <- apply(rbind(labs[allcombs[2,]],
	labs[allcombs[1,]]), 2,
	function(x)paste(x, collapse=" - "))

names(A.list) <- list.labs
dat <- data.frame(
	fit = do.call("c", effs),
	se.fit = do.call("c", se.effs),
	x = rep(quantseq, length(A.list)),
	contrast = rep(names(A.list), each=n), 
	stringsAsFactors=TRUE)
level <- level + ((1-level)/2)
dat$lower <- dat$fit - qt(level,
	obj$df.residual)*dat$se.fit
dat$upper <- dat$fit + qt(level,
	obj$df.residual)*dat$se.fit
res <- dat
for(i in c(1,2,3,5,6)){
    res[,i] <- round(res[,i], digits)
}
if(onlySig){
	sigs <- do.call(rbind,
		by(dat[,c("lower", "upper")],
		list(dat$contrast), function(x)
		c(max(x[,1]), min(x[,2]))))
	notsig <- which(sigs[,1] < 0 & sigs[,2] > 0)
    res <- res[-which(dat$contrast %in% names(notsig)), ]
}

if(type == "facs"){
	if(!plot){
        res <- list(out = res, mainvar = facvar, givenvar = quantvar)
        class(res) <- "iqq"
        print(res)
	}
	if(plot){
		rl <- range(c(res[, c("lower", "upper")]))
		if(rug)rl[1] <- rl[1] - (.05*length(faclevs))*diff(rl)
		p <- xyplot(fit ~ x | contrast, data=res, xlab = quantvar, ylab = "Predicted Difference", ylim = rl,
			lower=res$lower, upper=res$upper,
			prepanel = prepanel.ci, zl=TRUE,
		panel=function(x,y,subscripts,lower,upper,zl){
			panel.lines(x,y,col="black")
			if(ci){
				panel.lines(x,lower[subscripts], col="black", lty=2)
				panel.lines(x,upper[subscripts], col="black", lty=2)
			}
			if(zl)panel.abline(h=0, lty=3, col="gray50")
			if(rug){
				panel.doublerug(xa=l[[packet.number()]][[1]],xb=l[[packet.number()]][[2]])
		}
}
)
#	plot(p)
	return(p)
}
}
if(type == "slopes"){
	if(!plot){
	gq1 <- grep(paste(".*\\:", quantvar, "$", sep=""), names(b))
	gq2 <- grep(paste("^", quantvar, ".*\\:", sep=""), names(b))
	{if(length(gq1) == 0){qint <- gq2}
	else{qint <-  gq1}}
	if(length(qint) == 0){stop("Problem finding interaction coefficients")}
	qint <- c(grep(paste("^", quantvar, "$", sep=""), names(b)), qint)


	W <- matrix(0, nrow=length(faclevs), ncol=length(b))
	W[, qint[1]] <- 1
	W[cbind(1:length(faclevs), qint)] <- 1

	V <- vcov(obj)
	qeff <- c(W %*% b)
	qvar <- W %*% V %*% t(W)
	qse <- c(sqrt(diag(qvar)))
	qtstats <- c(qeff/qse)
	qpv <- c(2*pt(abs(qtstats), obj$df.residual, lower.tail=FALSE))

	qres <- sapply(list(qeff, qse, qtstats, qpv), function(x)sprintf("%.3f", x))
	colnames(qres) <- c("B", "SE(B)", "t-stat", "Pr(>|t|)")
	rownames(qres) <- faclevs
	names(qeff) <- faclevs
	res <- list(out = data.frame(eff = qeff, se = qse, tstat=qtstats, pvalue=qpv, stringsAsFactors=TRUE), varcor = qvar, mainvar = quantvar, givenvar = facvar)
    class(res) <- "iqq"
	print(res)
}
if(plot){
	intterm <- NULL
	if(paste(facvar, quantvar, sep=":") %in% colnames(attr(terms(obj), "factors"))){
		intterm <- paste(facvar, quantvar, sep="*")
	}
	if(paste(quantvar, facvar, sep=":") %in% colnames(attr(terms(obj), "factors"))){
		intterm <- paste(quantvar, facvar, sep="*")
	}
	if(is.null(intterm)){
		stop("No interaction in model\n")
	}
	e <- do.call(effect, c(list(term=intterm, mod=obj, default.levels=n, ...)))
	le <- as.list(by(mf[[quantvar]], list(mf[[facvar]]), function(x)x))

	edf <- data.frame(fit = e$fit, x = e$x[,quantvar],
	   fac = e$x[,facvar], se = e$se, stringsAsFactors=TRUE)
	edf$lower <- edf$fit - qt(.975, obj$df.residual)*edf$se
	edf$upper <- edf$fit + qt(.975, obj$df.residual)*edf$se

	yl <- range(c(edf$upper, edf$lower))
	xl <- range(edf$x) + c(-1,1)*.01*diff(range(edf$x))
	if(rug)yl[1] <- yl[1] - (.05*length(faclevs))*diff(yl)

	p <- xyplot(fit ~ x, group = edf$fac, data=edf,
		lower=edf$lower, upper=edf$upper,
		ylim = yl, xlim=xl,
		xlab = quantvar, ylab="Predicted Values",
		key=simpleKey(faclevs, lines=TRUE, points=FALSE),
		panel = function(x,y,groups, lower, upper, ...){
		if(ci){
			panel.transci(x,y,groups,lower,upper, ...)
		}
		else{
			panel.superpose(x=x,y=y, ..., panel.groups="panel.xyplot", type="l", groups=groups)
		}
		if(rug){
			for(i in 1:length(faclevs)){
				st <- (0) + (i-1)*.03
				end <- st + .02
				panel.rug(x=le[[i]],y=NULL,col=trellis.par.get("superpose.line")$col[i], start=st, end=end)

				}
			}
		})
	# plot(p)
	return(p)
}
}
}



#' Lattice panel function for translucent confidence intervals
#' 
#' This panel function is defined to plot translucent confidence intervals in a
#' single-panel, grouped (i.e., superposed) lattice display. Note, both lower
#' and upper must be passed directly to \code{xyplot} as they will be passed
#' down to the panel function.
#' 
#' 
#' @param x,y Data from the call to \code{xyplot}.
#' @param groups Variable used to created the superposed panels.
#' @param lower,upper 95\% lower and upper bounds of \code{y}.
#' @param ca Value of the alpha channel in [0,1]
#' @param ... Other arguments to be passed down to the plotting functions.
#' 
#' @export
#' 
#' @author Dave Armstrong
panel.transci <- function(x,y,groups,lower,upper, ca=.25, ...){
	ungroup <- unique(groups)
	sup.poly <- trellis.par.get("superpose.polygon")
	sup.poly.rgb <- col2rgb(sup.poly$col)/255
	sup.poly.rgb <- rbind(sup.poly.rgb, ca)
	rownames(sup.poly.rgb)[4] <- "alpha"
	ap <- apply(sup.poly.rgb, 2, function(x)lapply(1:4, function(y)x[y]))
	sup.poly$col <- sapply(ap, function(x)do.call(rgb, x))
    sup.poly$col <- rep(sup.poly$col, ceiling(length(ungroup)/length(sup.poly$col)))
	sup.line <- trellis.par.get("superpose.line")
    sup.line$col <- rep(sup.line$col, ceiling(length(ungroup)/length(sup.line$col)))
    sup.line$lwd <- rep(sup.line$lwd, ceiling(length(ungroup)/length(sup.line$lwd)))
    sup.line$lty <- rep(sup.line$lty, ceiling(length(ungroup)/length(sup.line$lty)))
	for(i in 1:length(ungroup)){
	panel.polygon(
		x=c(x[groups == ungroup[i]], rev(x[groups == ungroup[i]])),
		y = c(lower[groups == ungroup[i]], rev(upper[groups == ungroup[i]])),
		col = sup.poly$col[i], border="transparent")
	panel.lines(x[groups == ungroup[i]], y[groups == ungroup[i]], col=sup.line$col[i], lwd=sup.line$lwd[i], lty= sup.line$lty[i])
	}
}





#' Lattice panel function for two rug plots
#' 
#' This panel function is defined to plot two rugs, one on top of the other in
#' a multi-panel lattice display.
#' 
#' 
#' @param xa,xb Numeric vectors to be plotted.
#' @param regular Logical flag indicating whether rug is to be drawn on the
#' usual side (bottom/left) as opposed to the other side (top/right).
#' @param start,end Start and end points for the rug ticks on the y-axis.
#' @param x.units Character vectors, replicated to be of length two. Specifies
#' the (grid) units associated with start and end above. x.units are for the
#' rug on the x-axis and y-axis respectively (and thus are associated with
#' start and end values on the y and x scales respectively). See
#' \code{panel.rug} for more details.
#' @param lty,lwd Line type and width arguments (see \code{par} for more
#' details).
#' 
#' @export
#' 
#' @author Dave Armstrong
panel.doublerug <- function (xa = NULL, xb = NULL,
	regular = TRUE, start = if (regular) 0 else 0.97,
    end = if (regular) 0.03 else 1, x.units = rep("npc", 2),
    lty = 1, lwd = 1)
{
    x.units <- rep(x.units, length.out = 2)
    grid.segments(x0 = unit(xa, "native"), x1 = unit(xa, "native"),
        y0 = unit(start, x.units[1]), y1 = unit(end*.75, x.units[2]),
		gp=gpar(col.line="black", lty=lty, lwd=lwd))
	grid.segments(x0 = unit(xb, "native"), x1 = unit(xb, "native"),
		y0 = unit(end*1.25, x.units[1]), y1 = unit((2*end), x.units[2]),
		gp=gpar(col.line="black", lty=lty, lwd=lwd))
}



#' Lattice panel function for confidence intervals
#' 
#' This panel function is defined to plot confidence intervals in a multi-panel
#' lattice display. Note, both lower and upper must be passed directly to
#' \code{xyplot} as they will be passed down to the prepanel function.
#' 
#' 
#' @param x,y Data from the call to \code{xyplot}.
#' @param subscripts Variable used to created the juxtaposed panels.
#' @param lower,upper 95\% lower and upper bounds of \code{y}.
#' @param zl Logical indicating whether or not a horizontal dotted line at zero
#' is desired.
#' 
#' @export
#' 
#' @author Dave Armstrong
panel.ci <- function(x,y,subscripts,lower,upper,zl){
	panel.lines(x,y,col="black")
	panel.lines(x,lower[subscripts], col="black", lty=2)
	panel.lines(x,upper[subscripts], col="black", lty=2)
	if(zl)panel.abline(h=0, lty=3, col="gray50")
}



#' Lattice prepanel function for confidence intervals
#' 
#' This prepanel function is defined so as to allow room for all confidence
#' intervals plotted in a lattice display. Note, both lower and upper must be
#' passed directly to \code{xyplot} as they will be passed down to the prepanel
#' function.
#' 
#' 
#' @param x,y Data from the call to \code{xyplot}.
#' @param subscripts Variable used to created the juxtaposed panels.
#' @param lower,upper 95\% lower and upper bounds of \code{y}.
#' @return A list giving the ranges and differences in ranges of x and the
#' lower and upper bounds of y.
#' 
#' @export
#' 
#' @author Dave Armstrong
prepanel.ci <- function(x,y,subscripts, lower,upper){
    x2 <- as.numeric(x)
    list(xlim = range(x2, finite = TRUE),
         ylim = range(c(lower[subscripts], upper[subscripts]),
            finite = TRUE),
         dx = diff(range(x2, finite = TRUE)),
         dy = diff(range(c(lower[subscripts], upper[subscripts]),
            finite = TRUE)))
}


#' Lattice panel function for confidence intervals with capped bars
#' 
#' This panel function is defined to plot confidence intervals in a multi-panel
#' lattice display where the x-variable is categorical. Note, both lower and
#' upper must be passed directly to \code{xyplot} as they will be passed down
#' to the panel function.
#' 
#' 
#' @param x,y Data from the call to \code{xyplot}.
#' @param subscripts Variable used to created the juxtaposed panels.
#' @param lower,upper 95\% lower and upper bounds of \code{y}.
#' @param length Length of the arrow head lines.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' library(lattice)
#' library(effects)
#' data(Duncan, package="carData")
#' Duncan$inc.cat <- cut(Duncan$income, 3)
#' mod <- lm(prestige~ inc.cat * type + education,
#'   data=Duncan)
#' e1 <- effect("inc.cat*type", mod)
#' update(plot(e1), panel=panel.2cat)
#' 
panel.2cat <- function(x,y,subscripts,lower,upper, length=.2){
	panel.points(x,y, pch=16, col="black")
	panel.arrows(x, lower[subscripts], x, upper[subscripts], code=3, angle=90, length=length)
}




#' Test of linearity for Component + Residual Plots
#' 
#' This function estimates a linear model and a loess model on the
#' component-plus-residual plot (i.e., a partial residual plot) for each
#' quantitative variable in the model.  The residual sums of squares for each
#' are used to calculate an F-test for each quantitative variable.
#' 
#' 
#' @param model A model object of class \code{lm}
#' @param adjust.method Adjustment method for multiple-testing procedure, using
#' \code{p.adjust} from \code{stats}.
#' @param cat Number of unique values below which numeric variables are
#' considered categorical for the purposes of the smooth.
#' @param var Character string indicating the term desired for testing.  If
#' left \code{NULL}, the default value, all numeric variables will be tested.
#' @param span.as Logical indicating whether the span should be automatically
#' selected through AICC or GCV
#' @param span Span to be passed down to the \code{loess} function if
#' \code{span.as=FALSE}.
#' @param ... Other arguments to be passed down to the call to \code{loess}.
#' @return A matrix with the following columns for each variable:
#' \item{RSSp}{Residual sum-of-squares for the parametric (linear) model.}
#' \item{RSSnp}{Residual sum-of-squares for the non-parametric (loess) model.}
#' \item{DFnum}{Numerator degrees of freedom for the F-test: tr(S)-(k+1).}
#' \item{DFdenom}{Denominator degrees of freedom for the F-test: n-tr(S)}
#' \item{F}{F-statistic} \item{p}{p-value, potentially adjusted for multiple
#' comparisons.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' mod <- lm(prestige ~ income + education + women, data=Prestige)
#' crTest(mod)
#' 
crTest <- function(model, adjust.method="none", cat = 5, var=NULL, span.as = TRUE, span = 0.75, ...){
    cl <- attr(terms(model), "dataClasses")
	cl <- cl[which(cl != "factor")]
    mf <- model.frame(model)
    tabs <- apply(mf, 2, table)
    lens <- sapply(tabs, length)
    cats <- names(mf)[which(lens <= cat)]
    if(length(intersect(names(cl), cats)) > 0){
        cl <- cl[-which(names(cl) %in% cats)]
    }
    terms <- predictor.names(model)
	terms <- intersect(terms, names(cl))
    if (any(attr(terms(model), "order") > 1)) {
        stop("C+R plots not available for models with interactions.")
    }
    if(!is.null(var)){
        terms <- intersect(terms, var)
    }
    if(length(terms) == 0){
        stop(paste0(var, " not in list of model terms"))
    }
	terms.list <- list()
	orders <- sapply(terms, function(x)df.terms(model, x))
	for(i in 1:length(terms)){
        tmp.x <- {
            if (df.terms(model, terms[i]) > 1) predict(model, type = "terms", term = terms[i])
            else model.matrix(model)[, terms[i]]
        }
    	if(!is.null(colnames(tmp.x))){colnames(tmp.x) <- "x"}
    	terms.list[[i]] <- data.frame(x=tmp.x,
    		y = residuals.glm(model, "partial")[,terms[i]], stringsAsFactors=TRUE)
    }
if(!span.as){
    lo.mods <- lapply(terms.list, function(z)loess(y ~ x, data=z, span=span, ...))
}
if(span.as){
    lo.mods <- lapply(terms.list, function(z)loess.as(z$x, z$y, ...))
}
lin.mods <- lapply(terms.list, function(z)lm(y ~ x, data=z))
n <- nrow(model.matrix(model))
lo.rss <- sapply(lo.mods, function(x)sum(residuals(x)^2))
lm.rss <- sapply(lin.mods, function(x)sum(residuals(x)^2))
d1a <- sapply(lo.mods, function(x)x$one.delta)
d2a <- sapply(lo.mods, function(x)x$two.delta)
denom.df <- d1a^2/d2a
num.df <- (n-denom.df) - sapply(lin.mods, function(x)x$rank)
F.stats <- ((lm.rss-lo.rss)/num.df)/(lo.rss/denom.df)
pvals <- p.adjust(pf(F.stats, num.df, denom.df, lower.tail=FALSE), method=adjust.method)
out <- data.frame(
	RSSp = sprintf("%.2f", lm.rss),
	RSSnp = sprintf("%.2f", lo.rss),
	DFnum = sprintf("%.3f", num.df),
	DFdenom = sprintf("%.3f", denom.df),
	F = sprintf("%.3f", F.stats),
	p = sprintf("%.3f", pvals), stringsAsFactors=TRUE)
rownames(out) <- terms
return(out)
}


#' Test of Span Parameter in linearity for Component + Residual Plots
#' 
#' This function performs \code{crTest} for a user-defined range of span
#' parameters, optionally allowing for multiple testing corrections in the
#' p-values.
#' 
#' 
#' @param model A model object of class \code{lm}
#' @param spfromto A vector of two values across which a range of \code{n} span
#' values will be generated and tested.
#' @param n Number of span parameters to test.
#' @param adjust.method Adjustment method for multiple-testing procedure, using
#' \code{p.adjust} from \code{stats}.
#' @param adjust.type String giving the values over which the multiple testing
#' correction will be performed.  Here, \sQuote{both} refers to a multiple
#' testing correction done over all span parameters and all variables in the
#' model.  \sQuote{within} means the multiple testing correction should be done
#' within each model, but not across the span parameters and \sQuote{across}
#' means that the multiple testing correction should be for each variable
#' across the various span parameters, but not across variables within the same
#' model.  \sQuote{none} refers to a pass-through option of no multiple testing
#' procedure.
#' @return A list with two elements: \item{x}{Sequence of span values used in
#' testing} \item{y}{p-values for each variable for each span parameter}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' mod <- lm(prestige ~ income + education + women, data=Prestige)
#' tmp <- crSpanTest(mod, c(.1, .9), adjust.method="holm", 
#' 	adjust.type="both")
#' matplot(tmp$x, tmp$y, type="l")
#' 
crSpanTest <-
function(model, spfromto, n=10, adjust.method = "none", adjust.type = c("none", "across", "within", "both")){
	span.seq <- seq(from=spfromto[1], to=spfromto[2],
		length=n)
	adjust.type <- match.arg(adjust.type)
	out.list <- list()
	for(i in 1:length(span.seq)){
		out.list[[i]] <- crTest(model, adjust.method="none",
			span = span.seq[i])
	}
	pvals <- sapply(out.list, function(x)
		as.numeric(as.character(x[,"p"])))
	if(!is.matrix(pvals)){
		pvals <- matrix(pvals, nrow=1)
	}
	if(adjust.type == "within"){
		pvals <- apply(pvals, 2, p.adjust,
			method=adjust.method)
	}
	if(adjust.type == "across"){
		pvals <- t(apply(pvals, 1, p.adjust,
			method=adjust.method))
	}
	if(adjust.type == "both"){
		pvals <- matrix(p.adjust(c(pvals),
			method=adjust.method),
			nrow=nrow(pvals))
	}
	rownames(pvals) <- predictor.names(model)
	ret = list(x=span.seq, y=t(pvals))
    return(ret)
}


#' Standardize quantitative variables in a data frame
#' 
#' This function standardizes quantitative variables in a data frame while
#' leaving the others untouched.  This leaves not only factors, but also binary
#' variables (those with values 0, 1, or NA).
#' 
#' 
#' @param data A data frame.
#' @return A data frame with standardized quantitative variables
#' 
#' @export
#' 
#' @author Dave Armstrong
scaleDataFrame <-
function(data){
classes <- sapply(1:ncol(data), function(x)class(data[,x]))
dummies <- apply(data, 2, function(x)prod(x %in% c(0,1,NA)))
nn.ind <- union(which(dummies == 1), which(classes == "factor"))
num.ind <- setdiff(1:ncol(data), nn.ind)
if(length(num.ind) == 0){stop("No Numeric Variables in Data Frame")}
num.dat <- data[,num.ind]
nonnum.dat <- data[,nn.ind]
if(is.null(dimnames(num.dat))){
	num.dat <- data.frame(num.dat, stringsAsFactors=TRUE)
	names(num.dat) <- names(data)[num.ind]
}
if(is.null(dimnames(nonnum.dat))){
	nonnum.dat <- data.frame(nonnum.dat, stringsAsFactors=TRUE)
	names(nonnum.dat) <- names(data)[nn.ind]
}
newdat <- cbind(as.data.frame(scale(num.dat)), as.data.frame(nonnum.dat))
newdat
}


#' Create LaTeX or CSV versions of an Object Produced by CrossTable
#' 
#' \code{outXT} takes the output from \code{CrossTable} in the \code{gmodels}
#' package and produces either LaTeX code or CSV file that can be imported into
#' word processing software.
#' 
#' 
#' @param obj A list returned by \code{CrossTable} from the \code{gmodels}
#' package.
#' @param count Logical indicating whether the cell frequencies should be
#' returned.
#' @param prop.r Logical indicating whether the row proportions should be
#' returned.
#' @param prop.c Logical indicating whether the column proportions should be
#' returned.
#' @param prop.t Logical indicating whether the cell proportions should be
#' returned.
#' @param col.marg Logical indicating whether the column marginals should be
#' printed.
#' @param row.marg Logical indicating whether the row marginals should be
#' printed.
#' @param digits Number of digits to use in printing the proportions.
#' @param type String where \code{word} indicates a CSV file will be produced
#' and \code{latex} indicates LaTeX code will be generated.
#' @param file Connection where the file will be written, if \code{NULL} the
#' output will only be written to the console
#' @return A file containing LaTeX Code or CSV data to make a table
#' 
#' @export
#' 
#' @author Dave Armstrong
outXT <- function(obj, count=TRUE, prop.r = TRUE, prop.c = TRUE, prop.t = TRUE,
	col.marg=TRUE, row.marg=TRUE, digits = 3, type = "word", file=NULL){
	if(!(type %in% c("word", "latex"))){stop("type must be one of 'word' or 'latex'")}
	tmp.list <- list()
	k <- 1
	if(count){tmp.list[[k]] <- matrix(sprintf("%.0f", c(t(obj$t))), ncol=nrow(obj$t)); k <- k+1}
	if(prop.r){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""),
		c(t(obj$prop.r))), ncol=nrow(obj$prop.r)); k <- k+1}
	if(prop.c){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""),
		c(t(obj$prop.c))), ncol=nrow(obj$prop.c)); k <- k+1}
	if(prop.t){tmp.list[[k]] <- matrix(sprintf(paste("%.", digits, "f", sep=""),
		c(t(obj$prop.t))), ncol=nrow(obj$prop.t)); k <- k+1}
	out <- do.call(rbind, tmp.list)
	mat <- matrix(c(out), ncol=nrow(tmp.list[[1]]), byrow=TRUE)
	rownames(mat) <- NULL
	rn <- rep("", length=nrow(mat))
	rn[seq(1,nrow(mat), by=length(tmp.list))] <- rownames(obj$t)
	mat <- cbind(rn, mat)
	colnames(mat) <- c("", colnames(obj$t))
	if(col.marg){mat <- rbind(mat, c("Total", as.character(apply(obj$t, 2, sum))))}
	if(row.marg){
		rmarg <- rep("", nrow(mat))
		rmarg[seq(1, length(rmarg)-as.numeric(col.marg), by=length(tmp.list))] <- as.character(apply(obj$t, 1, sum))
		mat <- cbind(mat, rmarg)
		colnames(mat)[ncol(mat)] <- "Total"
	}
	if(row.marg & col.marg){mat[nrow(mat), ncol(mat)] <- sum(c(obj$t))}
	sink(file=ifelse(is.null(file), "tmp_xt.txt", file), type="output", append=FALSE)
	print(xtable(mat), include.rownames=FALSE, include.colnames=TRUE)
	sink(file=NULL)
	rl <- readLines(ifelse(is.null(file), "tmp_xt.txt", file))
	if(type == "word"){
		rl <- rl[-c(1:6,8)]
		rl <- rl[-c((length(rl)-3):length(rl))]
		rl <- gsub("\\\\", "", rl, fixed=TRUE)
		rl <- gsub(" & ", ",", rl, fixed=TRUE)
		if(!is.null(file)){writeLines(rl, con=file)}
	}
	print(rl, quote=FALSE)
}



#' Maximal First Differences for Generalized Linear Models
#' 
#' For objects of class \code{glm}, it calculates the change in predicted
#' responses, for discrete changes in a covariate holding all other variables
#' at their observed values.
#' 
#' The function calculates the average change in predicted probabiliy for a
#' discrete change in a single covariate with all other variables at their
#' observed values, for objects of class \code{glm}.  This function works with
#' polynomials specified with the \code{poly} function.
#' 
#' @aliases glmChange2
#' @param obj A model object of class \code{glm}.
#' @param varname Character string giving the variable name for which average
#' effects are to be calculated.
#' @param data Data frame used to fit \code{object}.
#' @param change A string indicating the difference in predictor values to
#' calculate the discrete change.  \code{sd} gives plus and minus one-half
#' standard deviation change around the median and \code{unit} gives a plus and
#' minus one-half unit change around the median.
#' @param R Number of simulations to perform.
#' @return \item{res}{A vector of values giving the average and 95 percent
#' confidence bounds} \item{ames}{The average change in predicted probability
#' (across all N observations) for each of the R simulations.}
#' \item{avesamp}{The average change in predicted probability for each of the N
#' observation (across all of the R simulations). }
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2), data=france, family=binomial)
#' glmChange2(left.mod, "age", data=france, "sd")
#' 
glmChange2 <-
function (obj, varname, data, change=c("unit", "sd"), R=1500)
{
    vars <- names(attr(terms(obj), "dataClasses"))[-1]
    vars <- gsub("poly\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("bs\\((.*?),.*?\\)", "\\1", vars)
    vars <- gsub("log\\((.*?),.*?\\)", "\\1", vars)
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
    b <- mvrnorm(R, coef(obj), vcov(obj))
    change <- match.arg(change)
    if(!is.factor(data[[varname]])){
    delt <- switch(change, unit = 1, sd = sd(data[[varname]],
        na.rm = TRUE))
    }else{
      delt <- 1
    }

    if (is.numeric(data[[varname]])) {
        d0 <- d1 <- data
        tmp0 <- d0[[varname]] - (0.5 * delt)
        tmp0 <- ifelse(tmp0 < min(d0[[varname]], na.rm=TRUE), min(d0[[varname]], na.rm=TRUE), tmp0)
        tmp1 <- d1[[varname]] + (0.5 * delt)
        tmp1 <- ifelse(tmp1 > max(d1[[varname]], na.rm=TRUE), max(d1[[varname]], na.rm=TRUE), tmp1)
        d0[[varname]] <- tmp0
        d1[[varname]] <- tmp1
        X0 <- model.matrix(obj, data = d0)
        X1 <- model.matrix(obj, data = d1)
        p0 <- family(obj)$linkinv(X0 %*% t(b))
        p1 <- family(obj)$linkinv(X1 %*% t(b))
        diff <- p1 - p0
        eff <- colMeans(diff)
        res <- matrix(c(mean(eff), quantile(eff, c(0.025, 0.975))),
            nrow = 1)
        colnames(res) <- c("mean", "lower", "upper")
        rownames(res) <- varname
        outres = list(res = res, ames=eff, avesamp = rowMeans(diff))
    }
    if (!is.numeric(data[[varname]]) & length(unique(na.omit(data[[varname]]))) ==
        2) {
        l <- obj$xlevels[[varname]]
        D0 <- D1 <- data
        D0[[varname]] <- factor(1, levels=1:2, labels=l)
        D1[[varname]] <- factor(2, levels=1:2, labels=l)
        X0 <- model.matrix(formula(obj), data=D0)
        X1 <- model.matrix(formula(obj), data=D1)
        p0 <- family(obj)$linkinv(X0 %*% t(b))
        p1 <- family(obj)$linkinv(X1 %*% t(b))
        diff <- p1 - p0
        eff <- colMeans(diff)
        res <- matrix(c(mean(eff), quantile(eff, c(0.025, 0.975))),
            nrow = 1)
        colnames(res) <- c("mean", "lower", "upper")
        rownames(res) <- varname
        outres = list(res = res, ames=eff, avesamp = rowMeans(diff))

    }
    if (!is.numeric(data[[varname]]) & length(unique(na.omit(data[[varname]]))) >
        2) {
        l <- obj$xlevels[[varname]]
        X.list <- list()
        for (j in 1:length(l)) {
            tmp <- data
            tmp[[varname]] <- factor(j, levels=1:length(l), labels=l)
            X.list[[j]] <- model.matrix(formula(obj), data=tmp)
        }
        combs <- combn(length(X.list), 2)
        d.list <- list()
        for (j in 1:ncol(combs)) {
            d.list[[j]] <- family(obj)$linkinv(X.list[[combs[2,
                j]]] %*% t(b)) - family(obj)$linkinv(X.list[[combs[1,
                j]]] %*% t(b))
        }
        eff <- sapply(d.list, colMeans)
        res <- apply(eff, 2, function(x) c(mean(x), quantile(x,
            c(0.025, 0.975))))
        cl <- array(l[combs], dim = dim(combs))
        rownames(res) <- c("mean", "lower", "upper")
        colnames(res) <- apply(cl[c(2, 1), ], 2, paste, collapse = "-")
        res <- t(res)
        outres = list(res = res, ames=eff, avesamp = sapply(d.list, rowMeans))
    }
    class(outres) <- "glmc2"
    print(outres)
}


#' Average Effect Plot for Generalized Linear Models
#' 
#' For objects of class \code{glm}, it calculates the change the average
#' predicted probability (like the one calculated by \code{glmChange2}) for a
#' hypothetical candidate set of values of a covariate.
#' 
#' The function plots the average effect of a model covariate, for objects of
#' class \code{glm}.  The function does not work with \code{poly} unless the
#' coefficients are provided as arguments to the command in the model (see
#' example below).
#' 
#' @param obj A model object of class \code{glm}.
#' @param varname Character string giving the variable name for which average
#' effects are to be calculated.
#' @param data Data frame used to fit \code{object}.
#' @param R Number of simulations to perform.
#' @param nvals Number of evaluation points at which the average probability
#' will be calculated.
#' @param plot Logical indicating whether plot should be returned, or just data
#' (if \code{FALSE}).
#' @param returnSim Logical indicating whether simulated predicted
#' probabilities should be returned.
#' @param ... Other arguments to be passed down to \code{xyplot}.
#' @return A plot or a data frame
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' p <- poly(france$lrself, 2)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2, coefs=attr(p, "coefs")), data=france, family=binomial)
#' aveEffPlot(left.mod, "age", data=france, plot=FALSE)
#' 
aveEffPlot <- function (obj, varname, data, R=1500, nvals=25, plot=TRUE, returnSim=FALSE, ...)
{
    vars <- all.vars(formula(obj))[-1]
    if(any(!(vars %in% names(data)))){
        vars <- vars[-which(!vars %in% names(data))]
    }
    rn <- vars
    var.classes <- sapply(vars, function(x) class(data[[x]]))
	b <- mvrnorm(R, coef(obj), vcov(obj), empirical=TRUE)
	if(is.numeric(data[[varname]])){
		s <- seq(min(data[[varname]], na.rm=TRUE), max(data[[varname]], na.rm=TRUE), length=nvals)
		dat.list <- list()
		for(i in 1:length(s)){
			dat.list[[i]] <- data
			dat.list[[i]][[varname]] <- s[i]
		}
		mm <- lapply(dat.list, function(x)model.matrix(obj, data=x))
		probs <- lapply(mm, function(x)family(obj)$linkinv(x %*% t(b)))
		cmprobs <- sapply(probs, colMeans)
		ciprobs <- t(apply(cmprobs, 2, function(x)c(mean(x), quantile(x, c(.025,.975)))))
		colnames(ciprobs) <- c("mean", "lower", "upper")
		ciprobs <- cbind(s, ciprobs)
		tmp <- as.data.frame(ciprobs)
		if(plot){
			pl <- xyplot(mean ~ s, data=tmp,  xlab=varname, ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				panel = function(x,y, lower, upper){
					panel.lines(x,y, col="black", lty=1)
					panel.lines(x, lower, col="black", lty=2)
					panel.lines(x, upper, col="black", lty=2)
			})
			return(pl)
		}
		else{
            if(returnSim){
                colnames(cmprobs) <- s
                class(cmprobs) <- "sims"
                out <- list(cis = tmp, probs=cmprobs)
                return(out)
            }
            else{
		    	return(tmp)
            }
        }
	}
	if(!is.numeric(data[[varname]])){
		l <- obj$xlevels[[varname]]
		dat.list <- list()
		for(j in 1:length(l)){
            dat.list[[j]] <- data
            dat.list[[j]][[varname]] <- factor(rep(j, nrow(data)), levels=1:length(l), labels=l)
			dat.list[[j]] <- model.matrix(formula(obj), data=dat.list[[j]])
		}
		probs <- lapply(dat.list, function(x)family(obj)$linkinv(x %*% t(b)))
		cmprobs <- sapply(probs, colMeans)
		ciprobs <- t(apply(cmprobs, 2, function(x)c(mean(x), quantile(x, c(.025,.975)))))
		colnames(ciprobs) <- c("mean", "lower", "upper")
		tmp <- as.data.frame(ciprobs)
		tmp$s <- factor(1:length(l), labels=l)
		if(plot){
			pl <- xyplot(mean ~ s, data=tmp, xlab="", ylab="Predicted Value", ...,
				lower=tmp$lower, upper=tmp$upper,
				prepanel = prepanel.ci,
				scales=list(x=list(at=1:length(l), labels=l)),
				panel = function(x,y, lower, upper){
					panel.points(x,y, col="black", lty=1, pch=16)
					panel.segments(x, lower, x, upper, lty=1, col="black")
			})
			return(pl)
		}
		else{
            if(returnSim){
                colnames(cmprobs) <- l
                class(cmprobs) <- "sims"
                out <- list(cis = tmp, probs=cmprobs)
                return(out)
            }
            else{
		    	return(tmp)
            }
        }
	}
}



#' AIC and BIC selection of number of spline knots
#' 
#' Calculates AIC and BIC for the selection of knots in a spline over values
#' (potentially including polynomials) up to a user-defined maximum.
#' 
#' 
#' @param form A formula detailing the model for which smoothing is to be
#' evaluated.
#' @param var A character string identifying the variable for which smoothing
#' is to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param includePoly Include linear and polynomial models up to, and including
#' \code{degree}-th order polynomials.
#' @param plot Logical indicating whether a plot should be returned.
#' @param criterion Statistical criterion to minimize in order to find the best
#' number of knots - AIC, BIC or Cross-validation.
#' @param cvk Number of groups for cross-validation
#' @param cviter Number of iterations of cross-validation to average over.  10
#' is the default but in real-world applications, this should be somewhere
#' around 200.
#' @return A plot, if \code{plot=TRUE}, otherwise a data frame with the degrees
#' of freedom and corresponding fit measure.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' NKnots(prestige ~ education + type, var="income", data=na.omit(Prestige), plot=FALSE)
#' 
NKnots <- function(form, var, data, degree=3, min.knots=1,
   max.knots=10, includePoly = FALSE, plot=FALSE, criterion=c("AIC", "BIC", "CV"),
   cvk=10, cviter=10){
   crit <- match.arg(criterion)
   k <- seq(min.knots, max.knots, by=1)
   forms <- vector("list", ifelse(includePoly, length(k)+3, length(k)))
   m <- 1
   if(includePoly){
      forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
      as.character(form)[3], " + ", var, sep=""))
      forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
      as.character(form)[3], "+ poly(", var,  ", 2)", sep=""))
      forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
      as.character(form)[3], "+ poly(", var,  ", 3)", sep=""))
      m <- 4
   }
   for(i in 1:length(k)){
      forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
      as.character(form)[3], "+ bs(", var, ", df=", degree+k[i],
        ", Boundary.knots=c(", min(data[[var]], na.rm=TRUE),", ", max(data[[var]], na.rm=TRUE), "))", sep=""))
      m <- m+1
   }
   if(crit %in% c("AIC", "BIC")){
       mods <- lapply(forms, function(x)lm(x, data=data))
       stats <- sapply(mods, function(x)do.call(crit, list(object=x)))
   }
   if(crit == "CV"){
      tmp.stats <- NULL
      for(j in 1:cviter){
       mods <- list()
       for(i in 1:length(forms)){
           tmp <- glm(forms[[i]], data=data, family=gaussian)
           tmpdat <- data[rownames(model.frame(tmp)), ]
           mods[[i]] <- cv.glm(tmpdat, tmp, K=cvk)
       }
       tmp.stats <- rbind(tmp.stats, sapply(mods, function(x)x$delta[1]))
     }
     stats <- colMeans(tmp.stats)
   }
   if(plot){
      k <- k+3
      if(includePoly){k <- c(1:3, k)}
      plot(k, stats, type="o", pch=16, col="black", xlab="# Degrees of Freedom", ylab = crit)
      points(k[which.min(stats)], min(stats), pch=16, col="red")
   }
}



#' Test of functional form assumption using B-splines
#' 
#' Estimate hypothesis test of lower- and higher-order non-linear relationships
#' against an assumed target relationship.
#' 
#' 
#' @param form A formula detailing the model for which smoothing is to be
#' evaluated.
#' @param var A character string identifying the variable for which smoothing
#' is to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param targetdf The assumed degrees of freedom against which the tests will
#' be conducted.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param adjust Method by which p-values will be adjusted (see
#' \code{\link{p.adjust}})
#' @return A matrix with the following columns: \item{F}{F statistics of test
#' of candidate models against target model} \item{DF1}{Numerator DF from
#' F-test} \item{DF2}{Denominator DF from F-test} \item{p(F)}{p-value from the
#' F-test} \item{Clarke}{Test statistic from the Clarke test}
#' \item{Pr(Better)}{The Clarke statistic divided by the number of
#' observations} \item{p(Clarke)}{p-value from the Clarke test.  (T) means that
#' the significant p-value is in favor of the Target model and (C) means the
#' significant p-value is in favor of the candidate (alternative) model.}
#' \item{Delta_AIC}{AIC(candidate model) - AIC(target model)}
#' \item{Delta_AICc}{AICc(candidate model) - AICc(target model)}
#' \item{Delta_BIC}{BIC(candidate model) - BIC(target model)}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' NKnotsTest(prestige ~ education + type, var="income", data=na.omit(Prestige), targetdf=3)
#' 
NKnotsTest <- function(form, var, data, targetdf = 1, degree=3, min.knots=1,
   max.knots=10, adjust="none"){
   k <- seq(min.knots, max.knots, by=1)
   forms <- vector("list", length(k)+3)
   m <- 1
   forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
   as.character(form)[3], " + ", var, sep=""))
   forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
   as.character(form)[3], "+ poly(", var,  ", 2)", sep=""))
   forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
   as.character(form)[3], "+ poly(", var,  ", 3)", sep=""))
   m <- 4
   for(i in 1:length(k)){
      forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
      as.character(form)[3], "+ bs(", var, ", df=", degree+k[i], ")", sep=""))
      m <- m+1
   }
   mods <- lapply(forms, function(x)lm(x, data=data))
   mods.df <- c(1:3, k+3)

   target.mod <- mods[[which(mods.df == targetdf)]]
   cand.mods <- mods
   cand.mods[[which(mods.df == targetdf)]] <- NULL
   tests <- lapply(cand.mods, function(x)as.matrix(anova(target.mod, x)))
   tests2 <- lapply(cand.mods, function(x)clarke_test(target.mod, x))
   num.df <- sapply(tests, function(x)abs(diff(x[,1])))
   denom.df <- sapply(tests, function(x)min(x[,1]))
   Fstats <- sapply(tests, function(x)x[2,5])
   pval <- p.adjust(sapply(tests, function(x)x[2,6]), method=adjust)
   cstats <- sapply(tests2, function(x)x$stat)
   cprobs <- sapply(tests2, function(x)x$stat/x$nobs)
   cminstat <- sapply(tests2, function(x)min(x$stat, x$nobs - x$stat))
   cbetter <-
   cp <- 2 * pbinom(cminstat, sapply(tests2, function(x)x$nobs), 0.5)
   pref <- sapply(tests2, function(x)ifelse(x$stat > x$nobs - x$stat,  "(T)", "(C)"))
   pref <- ifelse(cp > .05, "", pref)
   delta.aic <- sapply(cand.mods, AIC) - AIC(target.mod)
   delta.aicc <- sapply(cand.mods, AICc) - AICc(target.mod)
   delta.bic <- sapply(cand.mods, BIC) - BIC(target.mod)
   res <- cbind(Fstats, num.df, denom.df, pval, cstats, cprobs, cp, delta.aic, delta.aicc, delta.bic)
   sigchar <- ifelse(res[,4] < .05, "*", " ")
   sigchar2 <- ifelse(res[,7] < .05, "*", " ")
   strres <- NULL
   digs <- c(3,0,0,3, 0, 3, 3, 3, 3, 3)
   for(i in 1:10){
      tmp <- sprintf(paste("%.", digs[i], "f", sep=""), res[,i])
      if(i == 1){
         tmp <- paste(tmp, sigchar, sep="")
      }
      if(i == 5){
         tmp <- paste(tmp, sigchar2, sep="")
      }
      if(i == 7){
          tmp <- paste(tmp, pref,  sep=" ")
      }

      strres <- cbind(strres,tmp )
   }
   colnames(strres) <- c("F", "DF1", "DF2", "p(F)", "Clarke", "Pr(Better)", "p(Clarke)", "Delta_AIC", "Delta_AICc", "Delta_BIC")
   rownames(strres) <- paste("DF=", targetdf, " vs. DF=", mods.df[-targetdf], sep="")
   if(targetdf > 1){
      below <- strres[1:(targetdf-1), , drop=F]
      above <- strres[targetdf:nrow(strres),, drop=F]
      strres <- rbind(below, rep("", 10), above)
      rownames(strres)[targetdf] <- "   Target"
   }
   print(strres, quote=FALSE)
}



#' Significance Test for Loess vs. LM
#' 
#' Calculates an F test to evaluate significant differences between a LOESS
#' model and a parametric alternative estimated with \code{lm}
#' 
#' 
#' @param lmobj An object of class \code{lm}.
#' @param loessobj An object of class \code{loess}.
#' @param alpha Desired Type I error rate of test.
#' @return Printed output describing the results of the test.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' linmod <- lm(prestige ~ income, data=Prestige)
#' lomod <- loess(prestige ~ income, data=Prestige)
#' testLoess(linmod, lomod)
#' 
testLoess <- function(lmobj, loessobj, alpha=.05){
   n <- nrow(model.matrix(lmobj))
   if(n != loessobj$n){
       stop("Models estimated on different numbers of observations")
   }
rss0 <- sum(lmobj$residuals^2)
rss1 <- sum(loessobj$residuals^2)
d1a <- loessobj$one.delta; d2a <- loessobj$two.delta
dfdenom <- d1a^2/d2a
dfnum <- (n - dfdenom) - lmobj$rank
F0 <- ((rss0-rss1)/dfnum)/
    (rss1 / dfdenom)
cat("F = ", round(F0, 2), "\n", sep="")
pval <- pf(F0, dfnum, dfdenom, lower.tail=FALSE)
cat("Pr( > F) = ", round(pval, 2), "\n", sep="")
if(pval < alpha){
    cat("LOESS preferred to alternative\n")
}
else{
    cat("LOESS not statistically better than alternative\n")
}
}



#' Inspect a Variable in a Data Frame
#' 
#' Shows the variable label, factor levels (i.e., value labels) and frequency
#' distribution for the specified variable.
#' 
#' 
#' @aliases inspect inspect.data.frame inspect.tbl_df
#' @param data A data frame of class \code{data.frame} or \code{tbl_df}.
#' @param x A string identifying the name of the variable to be inspected.
#' @param ... Other arguments to be passed down, currently unimplemented.
#' @return A list with a variable label (if present), factor levels/value
#' labels (if present) and a frequency distribution
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' inspect(france, "vote")
#' 
inspect <- function(data, x, ...)UseMethod("inspect")

#' @method inspect tbl_df
#' @export
inspect.tbl_df <- function(data, x, ...){
  tmp <- data[[as.character(x)]]
  var.lab <- attr(tmp, "label")
  if(is.null(var.lab)){var.lab <- "No Label Found"}
  val.labs <- attr(tmp, "labels")
  if(is.null(val.labs)){val.labs <- sort(unique(tmp))}
  tab <- cbind(freq = table(tmp), prop = round(table(tmp)/sum(table(tmp), na.rm=TRUE), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}

##' @method inspect data.frame
##' @export
inspect.data.frame <- function(data, x, ...){
  var.lab <- attr(data, "var.label")[which(names(data) == x)]
  if(is.null(var.lab)){var.lab <- "No Label Found"}
  val.labs <- if(!is.null(levels(data[[x]]))){levels(data[[x]])}
    else {sort(unique(data[[x]]))}
  tab <- cbind(freq = table(data[[x]]), prop = round(table(data[[x]])/sum(table(data[[x]]), na.rm=TRUE), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}



#' Alternating Least Squares Optimal Scaling
#' 
#' Estimates the Alternating Least Squares Optimal Scaling (ALSOS) solution for
#' qualitative dependent variables.
#' 
#' \code{alsosDV} estimates the Alternating Least Squares Optimal Scaling
#' solution on the dependent variable.
#' 
#' @param form A formula with a dependent variable that will be optimally
#' scaled
#' @param data A data frame.
#' @param maxit Maximum number of iterations of the optimal scaling algorithm.
#' @param level Measurement level of the dependent variable 1=Nominal,
#' 2=Ordinal
#' @param process Nature of the measurement process: 1=discrete, 2=continuous.
#' Basically identifies whether tied observations will continue to be tied in
#' the optimally scaled variale (1) or whether the algorithm can untie the
#' points (2) subject to the overall measurement constraints in the model.
#' @param starts Optional starting values for the optimal scaling algorithm.
#' @param ... Other arguments to be passed down to \code{lm}.
#' @return A list with the following elements:
#' 
#' \item{result}{The result of the optimal scaling process}
#' 
#' \item{data}{The original data frame with additional columns adding the
#' optimally scaled DV}
#' 
#' \item{iterations}{The iteration history of the algorithm}
#' 
#' \item{form}{Original formula}
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
alsosDV <- function(form, data, maxit=30, level=2, process=1, starts=NULL,...){
	# process: 1 = discrete, 2 = continuous
	# level: 1 = nominal, 2 = ordinal
	# maxit: maximum number of optimal scaling iterations
	# source("http://www.quantoid.net/reg3/opscale_v14.R")

	rownames(data) <- 1:nrow(data)
	previous.rsquared <- niter <- 0
	rsquared.differ <- 1.0
	record <- c()
	form <- as.formula(form)
	mf <- model.frame(form, data)
	tmpdata <- data[which(rownames(data) %in% rownames(mf)),]
	dv <- form[[2]]
	tmpdata$dvar.os <- model.response(mf)
	if(!is.null(starts)){
		tmpdata$dvar.os <- starts
	}
	tmpmod <- lm(form, tmpdata, ...)
	while (rsquared.differ > .001 && niter <= maxit) {
	niter <- niter + 1
	reg.os <- update(tmpmod, dvar.os ~ .)
	rsquared.differ <- summary(reg.os)$r.squared - previous.rsquared
	previous.rsquared <- summary(reg.os)$r.squared
	record <- c(record, niter, summary(reg.os)$r.squared, rsquared.differ)
	   if (rsquared.differ > .001) {
	   dvar.pred <- predict(reg.os)
	   opscaled.dvar <- opscale(model.response(model.frame(form, data)), dvar.pred, level = level, process = process)
	   tmpdata$dvar.os <- opscaled.dvar$os
	   }
	}
	record <- matrix(round(record,4), ncol=3, byrow=TRUE)
	rownames(record) <- record[,1]
	record <- record[,-1]
	colnames(record) <- c("r-squared", "r-squared dif")
	tmpdata[[paste(dv, "_os", sep="")]] <- NA
	tmpdata[[paste(dv, "_os", sep="")]][match(rownames(tmpdata), rownames(mf))] <- tmpdata$dvar.os
	invisible(list(result=opscaled.dvar, data=tmpdata, iterations=record, formula=form))
}



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
#' @param alg Algorithm used to do sampling.  See \code{\link{stan}} for more
#' details.
#' @param ... Other arguments to be passed down to \code{stanfit}.
#' @return A list with the following elements:
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
balsos <- function(formula, data, iter=2500, chains = 1, alg = c("NUTS", "HMC", "Fixed_param"), ...){
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
fit <- stan(model_code=stancode, data = balsos.dat,
            iter = iter, chains = chains, algorithm=alg, ...)
ret <- list(fit = fit, y=y, X=X, form=formula)
class(ret) <- "balsos"
return(ret)
}

#figure out whether we need the col.alpha parameter in the panel.transci function.


#' Plot LOESS curve.
#' 
#' Plots the loess curve of the fitted values against a focal x-variable. .
#' 
#' Plots the fitted loess curve potentially with point-wise confidence bounds.
#' 
#' @param x An object of class \code{loess}.
#' @param ci Logical indicating whether point-wise confidence intervals should
#' be included around the fitted curve.
#' @param level The confidence level of the confidence intervals
#' @param linear Logical indicating whether the OLS line should also be
#' included.
#' @param addPoints Logical indicating whether or not points should be added to
#' the figure idtntifying the postion of individual observations.
#' @param col.alpha Value for alpha channel of the RGB color palette.
#' @param ... Other arguments to be passed down to \code{xyplot}.
#' @return A plot.
#' @method plot loess
#' 
#' @export
#' 
#' @author Dave Armstrong
plot.loess <- function(x, ..., ci=TRUE, level=.95, linear=FALSE, addPoints=FALSE, col.alpha=.5){
    alpha <- (1-level)/2
    tmp <- data.frame(x=x$x, y=x$y, stringsAsFactors=TRUE)
    xn <- as.character(formula(x$call)[[3]])
    yn <- as.character(formula(x$call)[[2]])
    names(tmp) <- c(xn, yn)
#    plot(formula(x$call), data=tmp, ...)
    xs <- seq(min(x$x), max(x$x), length=100)
    tmp2 <- data.frame(x=xs, y=rep(0, length(xs)), stringsAsFactors=TRUE)
    names(tmp2) <- names(tmp)
    preds <- predict(x, newdata=tmp2, se=TRUE)
    crit <- abs(qnorm(alpha))
    plot.dat <- data.frame(fit=preds$fit,
        lower=preds$fit - crit*preds$se.fit,
        upper = preds$fit + crit*preds$se.fit,
        x = xs,
        g = "loess", stringsAsFactors=TRUE)
    if(!linear){
     p <-    xyplot(fit ~ x, data=plot.dat, groups=plot.dat$g, ...,
            lower=plot.dat$lower,
            upper=plot.dat$upper,
            panel=panel.transci,
            prepanel=prepanel.ci)
    }
    if(linear){
        mod <- lm(formula(x$call), data=tmp)
        linpred <- predict(mod, newdata=tmp2, interval="confidence", level=level)
        linpred <- as.data.frame(linpred)
        names(linpred) <- c("fit", "lower", "upper")
        linpred$x <- xs
        linpred$g <- "linear"
        plot.dat <- rbind(plot.dat, linpred)
        p <- xyplot(fit ~ x, data=plot.dat, groups=plot.dat$g, ...,
            lower=plot.dat$lower,
            upper=plot.dat$upper,
            panel=panel.transci,
            prepanel=prepanel.ci)

    }
    print(p)
    if(addPoints){
        trellis.focus("panel",1,1)
        panel.points(x$x, x$y, col="black")
        trellis.unfocus()
    }
}


#' Estimate Derivatives of LOESS Curve.
#' 
#' Estimates the first derivatives of the LOESS curve.
#' 
#' 
#' @param obj An object of class \code{loess}.
#' @param delta Small change to be induced to estimate derivative.
#' @return A vector of first derivative values evaluated at each original
#' x-value.
#' 
#' @export
#' 
#' @author Dave Armstrong
loessDeriv <- function(obj, delta=.00001){
    newdf <- data.frame(y = obj$y, stringsAsFactors=TRUE)
    newdf[[obj$xnames[1]]] <- obj$x
    fit1 <- predict(obj, newdata=newdf, se=TRUE)
    newdf[[obj$xnames[1]]] <- obj$x + delta
    fit2 <- predict(obj, newdata=newdf, se=TRUE)
    deriv <- (fit2$fit - fit1$fit)/delta
    deriv
}



#' Regions of Statistical Significance in Interactions
#' 
#' Calculates the regions of statistical significance in interaction and
#' identifies the points at which the statistical significance of conditional
#' coefficients changes.
#' 
#' 
#' @param obj A model of class \code{glm} or class \code{lm}.
#' @param vars A character vector of the names of the two variables involved in
#' the interaction.
#' @param alpha Critical p-value of the test.
#' @return Printed output that identifies the change-points in statistical
#' significance.
#' 
#' @export
#' 
#' @author Dave Armstrong
changeSig <- function(obj, vars, alpha=.05){
    v1 <- which(names(obj$coef) == vars[1])
    v2 <- which(names(obj$coef) == vars[2])
    n12 <- ifelse(paste(vars, collapse=":") %in% names(coef(obj)), paste(vars, collapse=":"), paste(vars[2:1], collapse=":"))
    v12 <- which(names(obj$coef) == n12)

    B <- coef(obj)
    B.star <- B/qt(1-(alpha/2), obj$df.residual)
    V <- vcov(obj)

    # calculate solution for var1
    b <- 2*(V[v1,v12] - B.star[v1]*B.star[v12])
    a <- V[v12,v12] - B.star[v12]^2
    c <- V[v1,v1] - B.star[v1]^2

    x <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    x <- c(x, (-b - sqrt(b^2 - 4*a*c))/(2*a))
    x1 <- x
    est <- B[v1] + B[v12]*x
    v <- V[v1,v1] + x^2*V[v12,v12] + 2*x*V[v1,v12]
    s <- sqrt(v)
    lb1 <- est-qt(.975, obj$df.residual)*s
    ub1 <- est+qt(.975, obj$df.residual)*s


    # calculate solution for var2
    b <- 2*(V[v2,v12] - B.star[v2]*B.star[v12])
    a <- V[v12,v12] - B.star[v12]^2
    c <- V[v2,v2] - B.star[v2]^2

    x <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    x <- c(x, (-b - sqrt(b^2 - 4*a*c))/(2*a))
    x2 <- x
    est <- B[v2] + B[v12]*x
    v <- V[v2,v2] + x^2*V[v12,v12] + 2*x*V[v2,v12]
    s <- sqrt(v)
    lb2 <- est-qt(.975, obj$df.residual)*s
    ub2 <- est+qt(.975, obj$df.residual)*s

    tp1 <- mean(model.matrix(obj)[,v2] < x1[which.min(abs(lb1))])*100
    tp2 <- mean(model.matrix(obj)[,v2] < x1[which.min(abs(ub1))])*100
    tp3 <- mean(model.matrix(obj)[,v1] < x2[which.min(abs(lb2))])*100
    tp4 <- mean(model.matrix(obj)[,v1] < x2[which.min(abs(ub2))])*100
    lpc0 <- "(< Minimum Value in Data)"
    lpc100 <- "(> Maximum Value in Data)"
    lpcmid <- "(%.0fth pctile)"

    if(tp1 == 0){p1l <- lpc0}
    if(tp1 == 100){p1l <- lpc100}
    if(tp1 > 0 & tp1 < 100){p1l <- sprintf(lpcmid, tp1)}
    if(tp2 == 0){p1u <- lpc0}
    if(tp2 == 100){p1u <- lpc100}
    if(tp2 > 0 & tp2 < 100){p1u <- sprintf(lpcmid, tp2)}
    if(tp3 == 0){p2l <- lpc0}
    if(tp3 == 100){p2l <- lpc100}
    if(tp3 > 0 & tp3 < 100){p2l <- sprintf(lpcmid, tp3)}
    if(tp4 == 0){p2u <- lpc0}
    if(tp4 == 100){p2u <- lpc100}
    if(tp4 > 0 & tp4 < 100){p2u <- sprintf(lpcmid, tp4)}

    cat("LB for B(", vars[1], " | ", vars[2], ") = 0 when ", vars[2], "=", round(x1[which.min(abs(lb1))], 4), " ", p1l, "\n", sep="")
    cat("UB for B(", vars[1], " | ", vars[2], ") = 0 when ", vars[2], "=", round(x1[which.min(abs(ub1))], 4), " ", p1u, "\n", sep="")
    cat("LB for B(", vars[2], " | ", vars[1], ") = 0 when ", vars[1], "=", round(x2[which.min(abs(lb2))], 4), " ", p2l, "\n", sep="")
    cat("UB for B(", vars[2], " | ", vars[1], ") = 0 when ", vars[1], "=", round(x2[which.min(abs(ub2))], 4), " ", p2u, "\n", sep="")
    m1 <- round(rbind(x1, lb1, ub1), 5)
    colnames(m1) <- rep(vars[2], 2)
    m2 <- round(rbind(x2, lb2, ub2), 5)
    colnames(m2) <- rep(vars[1], 2)
    rownames(m1) <- rownames(m2) <- c("x", "Lower", "Upper")
    res <- list(m1, m2)
    names(res) <- vars
    invisible(res)
}



#' Testing Measurement Level Assumptions.
#' 
#' Uses the algorithm discussed in Armstrong and Jacoby 2018) to test the
#' intervality of a variable scaled by the Bayesian ALSOS method.
#' 
#' 
#' @param obj An object of class \code{balsos}.
#' @param cor.type Type of correlation to be used in the p-value calculation.
#' @return Printed output giving the Bayesian p-value evaluating the null that
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



#' Optimizing Yeo-Johnson Transformation
#' 
#' Uses \code{nlminb} to find the optimal Yeo-Johnson transformation parameters
#' conditional on a parametric model specification.
#' 
#' 
#' @param form A formula with a dependent variable that will be optimally
#' scaled
#' @param data A data frame.
#' @param trans.vars A character string identifying the variables that should
#' be transformed
#' @param round.digits Number of digits to round the transformation parameters.
#' @param ... Other arguments to be passed down to \code{lm}.
#' @return A linear model object that was estimated on the optimally
#' transformed variables.
#' 
#' @export
#' 
#' @author Dave Armstrong
yj_trans <- function(form, data, trans.vars, round.digits = 3, ...){
optim_yj <- function(pars, form, data, trans.vars, ...){
    form <- as.character(form)
    for(i in 1:length(trans.vars)){
        form <- gsub(trans.vars[i], paste("yeo.johnson(", trans.vars[i], ",", pars[i], ")", sep=""), form)
    }
    form <- as.formula(paste0(form[2], form[1], form[3]))
    m <- lm(form, data=data)
    ll <- logLik(m)
    -ll
}
    opt.pars <- nlminb(rep(1, length(trans.vars)), optim_yj, form=form, data=data, trans.vars = trans.vars, lower=0, upper=2)
    if(!is.null(round.digits)){
        opt.pars$par <- round(opt.pars$par, round.digits)
    }
    form <- as.character(form)
    for(i in 1:length(trans.vars)){
        form <- gsub(trans.vars[i], paste("yeo.johnson(", trans.vars[i], ",", opt.pars$par[i], ")", sep=""), form)
    }
    form <- as.formula(paste0(form[2], form[1], form[3]))
    out <- lm(form, data=data, ...)
    return(out)
}



#' Cross-validating Loess curve
#' 
#' Function provides the cross-validation error for the loess curve in a manner
#' that is amenable to optimization of the span.
#' 
#' 
#' @param span The span of the loess smoother.
#' @param form The formula that identifies the model
#' @param data A data frame containing the required variables.
#' @param cost Cost function to be passed down to loess.
#' @param K Number of folds for the cross-validation
#' @param numiter Number of times over which the cv error will be aggregated
#' @param which Return raw or corrected cv error
#' @return The cross-validation error from the loess curve.
#' 
#' @export
#' 
#' @author Dave Armstrong
cv.lo2 <- function (span, form, data, cost = function(y, yhat) mean((y - yhat)^2, na.rm=TRUE),
    K = n, numiter = 100, which=c("corrected", "raw")) {
    sample0 <- function (x, ...)x[sample.int(length(x), ...)]
    w <- match.arg(which)
    n <- nrow(data)
    f <- ceiling(n/K)
    tmp.y <- model.response(model.frame(as.formula(form), data=data))
	out.cv <- out.cost0 <- rep(NA, numiter)
	for(l in 1:numiter){
        s <- sample0(rep(1:K, f), n)
    	samps <- by(1:n, s, function(x)x)
        n.s <- sapply(samps, length)
        ms <- max(s)
        mod <- loess(formula=as.formula(form), data=data, span=span)
        cost.0 <- cost(tmp.y, fitted(mod))
        CV <- 0
        for (i in 1:length(samps)) {
            j.out <- samps[[i]]
            j.in <- setdiff(1:n, samps[[i]])
            d.mod <- loess(as.formula(form), data=data[j.in, ], span=span)
            p.alpha <- n.s[i]/n
            cost.i <-  cost(tmp.y[j.out], predict(d.mod, data[j.out,
                , drop = FALSE]))
            CV <- CV + p.alpha * cost.i
            cost.0 <- cost.0 - p.alpha * cost(tmp.y, predict(d.mod, data))
        }
	out.cv[l] <- CV
	out.cost0[l] <- cost.0
	names(out.cv) <- names(out.cost0) <- NULL
}
    return(switch(w, raw=mean(out.cv), corrected=mean(out.cv+out.cost0)))
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
#' @return A lattice graph produce by a call to \code{xyplot}.
#' @method plot balsos
#' 
#' @export
#' 
#' @author Dave Armstrong
plot.balsos <- function(x, ..., freq=TRUE, offset=.1){
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
    os = s$statistics[,1],
    lower = s$quantiles[,1],
    upper = s$quantiles[,5],
    orig = sort(unique(x$y), stringsAsFactors=TRUE)
)
tp <- trellis.par.get()$superpose.symbol
if(!freq){
xyplot(os ~ orig, data=newdf,
    panel = function(x,y, ...){
        panel.points(x,y, cex=.6, col=tp$col[1], pch=tp$pch[1])
        panel.segments(x, newdf$lower, x, newdf$upper, col="black", cols=tp$col[1])
    })
}
else{
    newdf$freq <- a1$result$os[mins[,2]]
    xyplot(os ~ orig, data=newdf,
        panel = function(x,y, ...){
            panel.points(x-offset,y, cex=.6, col=tp$col[1], pch=tp$pch[1])
            panel.segments(x-offset, newdf$lower, x-offset, newdf$upper, col=tp$col[1])
            panel.points(x+offset, newdf$freq, col=tp$col[2], pch=tp$pch[2])
        })
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

central <- function(x){
    if(is.factor(x)){
        tab <- table(x)
        m <- which.max(tab)
        cent <- factor(m, levels=1:length(levels(x)), labels=levels(x))
    }
    if(!is.factor(x) & length(unique(na.omit(x))) <= 10){
        tab <- table(x)
        cent <- as.numeric(names(tab)[which.max(tab)])
    }
    if(!is.factor(x) & length(unique(na.omit(x))) > 10){
        cent <- median(x, na.rm=TRUE)
    }
    return(cent)
}



#' Print Confidence Intervals for Predicted Probabilities and Their Differences
#' 
#' Print method for output from the \code{probci} function.
#' 
#' 
#' @param x A object of class \code{diffci} produced by \code{\link{probci}}.
#' @param digits How many digits to round output.
#' @param filter A named list of values where the names indicate the variable
#' to be filtered and the values in the vector indicate the values to include
#' for the filtering variable.
#' @param const A string identifying the name of the variable to be held
#' constant across comparisons.
#' @param onlySig Logical indicating whether all differes should be displayed
#' or only those significant at the 95\% two-tailed level.
#' @param ... Other arguments to be passed down to print, currently
#' unimplemented.
#' @return An data frame with the following variables: \item{variables}{The
#' variables and the values at which they are held constant.  For example,
#' \code{tmp1} would be the first value of \code{tmp} used in the probability
#' calculation and \code{tmp2} would be the second value of \code{tmp} used in
#' the probability calculation. } \item{pred_prob}{The difference in predicted
#' probability given the following change in \code{X}:
#' \code{tmp2}-\code{tmp1}.} \item{lower, upper}{The lower and upper 95\%
#' confidence bounds.}
#' @method print diffci
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2, raw=TRUE), data=france, family=binomial)
#' out <- probci(left.mod, france, numQuantVals=3, 
#'     changeX=c("retnat", "lrself"))
#' print(out, filter=list(retnat=c("Better", "Worse")))
#' print(out, filter=list(retnat=c("Better", "Worse")),
#'      const="lrself")
#' print(out, const="retnat")
#' 
print.diffci <- function(x, ..., digits=4, filter=NULL, const = NULL, onlySig=FALSE){
    D <- x[[2]]
    if(!is.null(filter)){
        vn <- names(filter)
        w <- array(dim=c(nrow(D), length(vn)))
        for(i in 1:length(filter)){
            w[,i] <- apply(D[,grep(paste("^", vn[i], "[1-2]", sep=""), names(D))], 1,
                function(x)(x[1] %in% filter[[i]] & x[2] %in% filter[[i]]))
        }
        w <- which(apply(w, 1, prod) == 1)
        if(length(w) == 0){
            stop("No observations matching filter conditions")
        }
        else{
            D <- D[w, ]
        }
    }
    if(onlySig){
        sig <- which(sign(D$lower) == sign(D$upper))
        if(length(sig) == 0){
            stop("No observations matching filter conditions that are also significant")
        }
        else{
            D <- D[sig, ]
        }
      }
    if(!is.null(const)){
      D <- D[which(D[[paste0(const, "1")]] == D[[paste0(const, "2")]]), ]
    }
    cnd <- colnames(D)
    D2 <- array(dim=dim(D))
    fmts <- c(paste0("%.", digits, "f"), "%.1i")
    for(i in 1:ncol(D)){
        if(is.numeric(D[[i]])){
            D2[,i] <- sprintf(ifelse(is.integer(D[[i]]) |
                all(round(D[[i]], 5) == as.integer(D[[i]])),
                fmts[2], fmts[1]), D[[i]])
        }
        if(is.factor(D[[i]]) | is.character(D[[i]])){
            D2[,i] <- as.character(D[[i]])
        }
    }
    colnames(D2) <- cnd
    rownames(D2) <- 1:nrow(D2)
    print.data.frame(as.data.frame(D2))
}




#' Confidence Intervals for Predicted Probabilities and Their Differences
#' 
#' Calculates predicted probabilities for any combination of x-variable values
#' holding all other variables constant at either typical values (average case
#' approach) or at observed values (average effect approach).
#' 
#' Calculates predicted probabilities for any combination of x-variable values
#' holding all other variables constant at either typical values (average case
#' approach) or at observed values (average effect approach).  The function
#' uses a parametric bootstrap to provide generate confidence bounds for
#' predicted probabilities and their differences.  The confidence intervals
#' produced are raw percentile interviews (at the 5\% level).
#' 
#' @param obj A model of class \code{glm}, particularly those with
#' \code{binomial} family.
#' @param data Data frame used to estimate \code{obj}.
#' @param .b A vector of coefficients to be passed down to the simulation.  If
#' \code{NULL}, \code{coef()} will be used to obtain coefficients from
#' \code{obj}.
#' @param .vcov A parameter variance covariance matrix to be passed to the
#' simulation.  If \code{NULL}, \code{vcov()} will be used to obtain the
#' variance-covariance matrix of the parameters.
#' @param changeX A vector of strings giving the names of variables for which
#' changes are desired.
#' @param numQuantVals For quantitative variables, if no x-values are specified
#' in \code{xvals}, then \code{numQuantVals} gives the number of values used
#' across the range of the variable.
#' @param xvals A named list of values used to make the predictions.  The names
#' in the list should correspond with the variable names specified in
#' \code{changeX}.
#' @param type Type of effect to be generated.  \code{aveEff} produces the
#' average first difference across all observed values of \code{X}, while
#' \code{aveCase} gives the first difference holding all other variables
#' constant at typical values.
#' @return An data frame with the following variables: \item{variables}{The
#' variables and the values at which they are held constant.  For example,
#' \code{tmp1} would be the first value of \code{tmp} used in the probability
#' calculation and \code{tmp2} would be the second value of \code{tmp} used in
#' the probability calculation. } \item{pred_prob}{The difference in predicted
#' probability given the following change in \code{X}:
#' \code{tmp2}-\code{tmp1}.} \item{lower, upper}{The lower and upper 95\%
#' confidence bounds.}
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' 
#' data(france)
#' left.mod <- glm(voteleft ~ male + age + retnat + 
#' 	poly(lrself, 2, raw=TRUE), data=france, family=binomial)
#' out <- probci(left.mod, france, changeX="retnat")
#' out
#' out2 <- probci(left.mod, france, changeX="lrself",
#'     xvals = list(lrself = c(1,10)))
#' out2
#' out3 <- probci(left.mod, france, changeX=c("lrself", "retnat"), 
#'     xvals = list(lrself = c(1,10)))
#' out3
#' 
probci <- function(obj, data, .b = NULL, .vcov=NULL, changeX=NULL, numQuantVals=5, xvals = NULL, type=c("aveEff", "aveCase")){
    type <- match.arg(type)
    vn <- changeX
    if(length(vn) == 0){stop("Need at least one variable to change")}
    vals <- vector(length=length(vn), mode="list")
    names(vals) <- vn
    for(i in 1:length(vn)){
        if(is.factor(data[[vn[i]]])){
            vals[[i]] <- factor(1:length(levels(data[[vn[i]]])), labels=levels(data[[vn[i]]]))
        }
        if(!is.null(xvals[[vn[i]]])){
          vals[[i]] <- xvals[[vn[i]]]
        }
        else{
        if(!is.factor(data[[vn[i]]]) & length(unique(na.omit(data[[vn[i]]]))) <= numQuantVals){
            vals[[i]] <- sort(unique(na.omit(data[[vn[i]]])))
        }
        if(!is.factor(data[[vn[i]]]) & length(unique(na.omit(data[[vn[i]]]))) > numQuantVals){
                vals[[i]] <- seq(min(data[[vn[i]]], na.rm=TRUE), max(data[[vn[i]]], na.rm=TRUE), length = numQuantVals)
        }
    }
}
    egvals <- do.call(expand.grid, vals)
    if(is.null(.b)){
        b <- coef(obj)
    }
    else{
        b <- .b
    }
    if(is.null(.vcov)){
        v <- vcov(obj)
    }
    else{
        v <- .vcov
    }
    bmat <- t(mvrnorm(2500, b, v, empirical=TRUE))
    if(type == "aveCase"){
        av <- all.vars(obj$formula)
        others <- av[-which(av %in% vn)]
        othervals <- lapply(1:length(others), function(i)central(data[[others[i]]]))
        names(othervals) <- others
        otherdat <- do.call(data.frame, c(othervals, stringsAsFactors=TRUE))
        alldat <- do.call(expand.grid, c(vals, otherdat))
        X <- model.matrix(formula(obj), data=alldat)
        probs <- t(family(obj)$linkinv(X %*% bmat))
        probci <- t(apply(probs, 2, quantile, c(.5,.025,.975)))[,,drop=FALSE]
        dc <- combn(nrow(X), 2)
        D <- matrix(0, nrow=nrow(X), ncol=ncol(dc))
        D[cbind(dc[1,], 1:ncol(dc))] <- -1
        D[cbind(dc[2,], 1:ncol(dc))] <- 1
        diffprobs <- probs %*% D
        ev <- sapply(1:ncol(egvals), function(i)paste(colnames(egvals)[i], "=", egvals[,i], sep=""))[,,drop=FALSE]
        probn <- apply(ev, 1, paste, collapse=", ")
        rownames(probci) <- probn
        ev2 <- egvals[dc[2,], , drop=FALSE]
        ev2 <- as.matrix(sapply(1:ncol(ev2), function(i)paste(colnames(ev2)[i], "=", ev2[,i], sep="")))
        n1 <- apply(ev2, 1, paste, collapse=", ")

        ev1 <- egvals[dc[1,], , drop=FALSE]
        ev1 <- as.matrix(sapply(1:ncol(ev1), function(i)paste(colnames(ev1)[i], "=", ev1[,i], sep="")))
        n2 <- apply(ev1, 1, paste, collapse=", ")
        n <- paste("(", n1,  ") - (", n2, ")", sep="")
        diffci <- t(apply(diffprobs, 2, quantile, c(.5,.025,.975)))[,,drop=FALSE]
        rownames(diffci) <- n
        colnames(diffci) <- colnames(probci) <- c("pred_prob", "lower", "upper")
        tmp1 <- egvals[dc[1,], ]
        tmp2 <- egvals[dc[2,], ]
        names(tmp1) <- paste(names(tmp1), "1", sep="")
        names(tmp2) <- paste(names(tmp2), "2", sep="")
        diffci <- cbind(tmp1, tmp2, as.data.frame(diffci))[,,drop=FALSE]
        probci <- as.data.frame(probci)
    }
    if(type == "aveEff"){
        probs <- vector(mode="list", length=nrow(egvals))
        for(i in 1:nrow(egvals)){
            tmp <- data
            for(j in 1:ncol(egvals)){
                tmp[, names(egvals)[j]] <- egvals[i,j]
            }
            X <- model.matrix(formula(obj), data=tmp)
            probs[[i]] <- family(obj)$linkinv(X %*% bmat)
        }
        probci <- t(apply(sapply(probs, colMeans), 2, quantile, c(.5,.025, .975)))[,,drop=FALSE]
        ev <- sapply(1:ncol(egvals), function(i)paste(colnames(egvals)[i], "=", egvals[,i], sep=""))[,,drop=FALSE]
        probn <- apply(ev, 1, paste, collapse=", ")
        rownames(probci) <- probn
        dc <- combn(nrow(egvals), 2)
        diffprobs <- matrix(NA, nrow=2500, ncol=ncol(dc))
        for(i in 1:ncol(dc)){
            diffprobs[,i] <- colMeans(probs[[dc[2,i]]] - probs[[dc[1,i]]])
        }
        ev2 <- egvals[dc[2,], , drop=FALSE]
        ev2 <- as.matrix(sapply(1:ncol(ev2), function(i)paste(colnames(ev2)[i], "=", ev2[,i], sep="")))
        n1 <- apply(ev2, 1, paste, collapse=", ")

        ev1 <- egvals[dc[1,], , drop=FALSE]
        ev1 <- as.matrix(sapply(1:ncol(ev1), function(i)paste(colnames(ev1)[i], "=", ev1[,i], sep="")))
        n2 <- apply(ev1, 1, paste, collapse=", ")
        n <- paste("(", n1,  ") - (", n2, ")", sep="")
        diffci <- t(apply(diffprobs, 2, quantile, c(.5,.025,.975)))[,,drop=FALSE]
        rownames(diffci) <- n
        colnames(diffci) <- colnames(probci) <- c("pred_prob", "lower", "upper")
        tmp1 <- egvals[dc[1,], ]
        tmp2 <- egvals[dc[2,], ]
        names(tmp1) <- paste(names(tmp1), "1", sep="")
        names(tmp2) <- paste(names(tmp2), "2", sep="")
        diffci <- cbind(tmp1, tmp2, as.data.frame(diffci))[,,drop=FALSE]
        probci <- as.data.frame(probci)
    }
    g <- grep("^tmp[1-2]", names(diffci))
    if(length(g) == 2 & length(changeX) == 1){
        names(diffci) <- gsub("tmp", changeX, names(diffci))
    }
    rownames(diffci) <- NULL
    res <- list("Predicted Probabilities"=probci, "Difference in Predicted Probabilities"=diffci, plot.data = cbind(egvals, probci))
    class(res) <- "diffci"
    return(res)
}


##' @method plot diffci
plot.diffci <- function(x, ..., xvar = NULL, condvars = NULL){
    D <- x$plot.data
    rg <- range(D[,c("lower", "upper")])
    yl <- rg +abs(diff(rg))*.1 %o% c(-1,1)
    Dnp <- names(D)[which(!(names(D) %in% c("pred_prob", "lower", "upper")))]
    if(is.null(xvar)){
        if(length(Dnp) == 1){
            xvar <- Dnp
        }
        if(length(Dnp) > 1){
            lens <- apply(D[,Dnp, drop=F], 2, function(x)length(unique(x)))
            xvar <- Dnp[which.max(lens)]
        }
    }
    if(is.null(condvars)){
        if(length(Dnp) > 1){
            condvars <- Dnp[-which(Dnp == xvar)]
        }
    }
    if(!is.null(condvars)){
        for(i in 1:length(condvars)){
            levels(D[[condvars[i]]]) <- paste(condvars[i], levels(D[[condvars[i]]]), sep=": ")
        }
        ncondvars <- length(condvars)
        condvars <- paste(condvars, collapse="+")
    }
    l.arrow <- 1/(2*length(unique(D[[xvar]])))
    form <- as.formula(paste("pred_prob ~ ", xvar,  ifelse(is.null(condvars), "", paste("|", condvars))))
    plotargs <- list(x = form, data=D, lower=D$lower, upper=D$upper,
        par.settings = list(strip.background = list(col="gray80")))
    pa2 <- list(...)
    if("length" %in% names(pa2)){
        l.arrow <- pa2[["length"]]
        pa2[[which(names(pa2) == "length")]] <- NULL
    }
    plotargs <- c(plotargs, pa2)
    if(!("zl" %in% names(plotargs))){
        plotargs$zl <- TRUE
    }
    if(!("ylim" %in% names(plotargs))){
        plotargs$ylim <- c(yl)
    }
    if(!("panel" %in% names(plotargs))){
        plotargs$panel <- function (x, y, subscripts, lower, upper, length=l.arrow){
        panel.points(x, y, pch = 16, col = "black")
        panel.arrows(x, lower[subscripts], x, upper[subscripts],
            code = 3, angle = 90, length = length)
}
    }

    p <- do.call(xyplot, plotargs)
    if(!is.null(condvars)){
        if(ncondvars >= 2){
            return(useOuterStrips(p))
        }
     }
     else{
         return(p)
     }
}



#' Simulated F-test for Linear Interactions
#' 
#' Simluates the sampling distribution of the F statistic when comparing a
#' linear intreraction model to a generalized additive model with a smooth over
#' the two variables in the interaction.
#' 
#' In simple simulations, an F-test of a linear interaction relative to a
#' smooth interaction with a GAM using a nominal .05 type I error rate, has an
#' actual type I error rate of more than double the nominal rate (this tended
#' to be in the low teens).  This function tries to build the F-distribution
#' using simulation.  First, it uses the coefficients from the linear
#' interaction model, multiplies them by the coefficients from the linear
#' interaction model and for each iteration of the simulation, it creates the
#' simulated dependent variable by adding a random error to the linear
#' predictor with the same standard deviation as the residual standard
#' deviation from the linear interaction model.  All of that is to say that
#' this model has all of the same features as the linear interaction model,
#' except that we are certain that this is the right model.  The algorithm then
#' estimates both the linear interaction model and the GAM with a smooth
#' interaction on the original X variables and the new simulated y variable.
#' The F-test is performed and the F-statistic saved for each iteraction.  The
#' algorithm then calculates the probability of being to the right of the
#' observed F-statistic in the simulated F-distribution.
#' 
#' @param m1 An object of class \code{gam} estimated with the \code{mgcv}
#' package.  This model should be linear in the interaction of the two
#' x-variables of interest.
#' @param m2 An object of class \code{gam} esimtated with the \code{mgcv}
#' package.  This model should contain a smooth interaction.  For two
#' continuous variables, this should be done with \code{te()} unless the
#' variables are measured in the same units (e.g., spatial coordinates) in
#' which case the usual thin-plate regression spline will work.  For
#' categorical moderators, you should use the \code{s(x, by=D0)} and \code{s(x,
#' by=D)} (for a dummy variable moderator, \code{D}, where \code{D0=1} when
#' \code{D=0}.  Remember to include \code{D} as a parametric term in the model
#' as well to account for the intercept difference between the two smooth
#' terms.)
#' @param data Data frame used to estimate both models
#' @param R Number of simulated F values to create.
#' @param ranCoef Logcial indicating whether the coefficients should be treated
#' as fixed or whether they should be drawn from their implied sampling
#' distribution for each iteration of the simulation.
#' @return \item{obsF}{The observed F-statistic from the test on the original
#' models.} \item{Fdist}{The \code{R} different F-statistics calculated at each
#' iteration of the simulation.}
#' 
#' @export
#' 
#' @author Dave Armstrong
testGAMint <- function(m1, m2, data, R=1000, ranCoef=FALSE){
    v1 <- all.vars(formula(m1))
    v2 <- all.vars(formula(m2))
    av <- union(v1, v2)
    tmp <- data[,av]
    tmp <- na.omit(tmp)
    b1 <- coef(m1)
    X <- model.matrix(m1)
    obsF <- anova(m1, m2, test="F")[2,5]
    Fstat <- rep(NA, R)
    i <- 1
    while(i <= R){
        if(ranCoef){
            b1 <- mvrnorm(1, coef(m1), vcov(m1))
        }
        if(all("lm" %in% class(m1))){
            sig <- summary(m1)$sigma
        }
        if("gam" %in% class(m1)){
            sig <- sqrt(summary(m1)$scale)
        }
        e <- rnorm(nrow(X), 0, sig)
        tmp$ynew <- X %*% b1 + e
        m2a <- update(m2, ynew ~ ., data=tmp)
        m1a <- update(m1, ynew ~ ., data=tmp)
        tmpF <- anova(m1a, m2a, test="F")[2,5]
        if(!is.na(tmpF)){
            Fstat[i] <- tmpF
            i <- i+1
        }
    }
    return(list(obsF = obsF, Fdist = Fstat))
}



#' Conditional Effects Plots for Interactions in Linear Models
#' 
#' Generates two conditional effects plots for two interacted continuous
#' covariates in linear models.
#' 
#' This function does the same thing as \code{\link{DAintfun2}}, but presents
#' effects only at the mean of the conditioning variable and the mean +/- 1
#' standard deviation.
#' 
#' @param obj A model object of class \code{lm}
#' @param varnames A two-element character vector where each element is the
#' name of a variable involved in a two-way interaction.
#' @param varcov A variance-covariance matrix with which to calculate the
#' conditional standard errors.  If \code{NULL}, it is calculated with
#' \code{vcov(obj)}.
#' @param name.stem A character string giving filename to which the appropriate
#' extension will be appended
#' @param xlab Optional vector of length two giving the x-labels for the two
#' plots that are generated.  The first element of the vector corresponds to
#' the figure plotting the conditional effect of the first variable in
#' \code{varnames} given the second and the second element of the vector
#' corresponds to the figure plotting the conditional effect of the second
#' variable in \code{varnames} conditional on the first.
#' @param ylab Optional vector of length two giving the y-labels for the two
#' plots that are generated.  The first element of the vector corresponds to
#' the figure plotting the conditional effect of the first variable in
#' \code{varnames} given the second and the second element of the vector
#' corresponds to the figure plotting the conditional effect of the second
#' variable in \code{varnames} conditional on the first.
#' @param plot.type One of \sQuote{pdf}, \sQuote{png}, \sQuote{eps} or
#' \sQuote{screen}, where the one of the first three will produce two graphs
#' starting with \code{name.stem} written to the appropriate file type and the
#' third will produce graphical output on the screen.
#' @return \item{graphs}{Either a single graph is printed on the screen (using
#' \code{par(mfrow=c(1,2))}) or two figures starting with \code{name.stem} are
#' produced where each gives the conditional effect of one variable based on
#' the values of another.}
#' @author Dave Armstrong
#' @references Brambor, T., W.R. Clark and M. Golder.  (2006) Understanding
#' Interaction Models: Improving Empirical Analyses.  Political Analysis 14,
#' 63-82.\cr Berry, W., M. Golder and D. Milton.  (2012) Improving Tests of
#' Theories Positing Interactions.  Journal of Politics.
#' 
#' @export
#' 
#' @examples
#' 
#' data(InteractionEx)
#' mod <- lm(y ~ x1*x2 + z, data=InteractionEx)
#' DAintfun3(mod, c("x1", "x2"))
#' 
DAintfun3 <-
function (obj, varnames, varcov=NULL, name.stem = "cond_eff",
	xlab = NULL, ylab=NULL, plot.type = "screen")
{
    rseq <- function(x) {
        rx <- range(x, na.rm = TRUE)
        seq(rx[1], rx[2], length = 25)
    }
    MM <- model.matrix(obj)
    v1 <- varnames[1]
    v2 <- varnames[2]
    ind1 <- grep(paste0("^",v1,"$"), names(obj$coef))
    ind2 <- grep(paste0("^",v2,"$"), names(obj$coef))
    indboth <- which(names(obj$coef) %in% c(paste0(v1,":",v2),paste0(v2,":",v1)))
    ind1 <- c(ind1, indboth)
    ind2 <- c(ind2, indboth)
    s1 <- c(mean(MM[,v1]) - sd(MM[,v1]), mean(MM[,v1]), mean(MM[,v1]) + sd(MM[,v1]))
    s2 <- c(mean(MM[,v2]) - sd(MM[,v2]), mean(MM[,v2]), mean(MM[,v2]) + sd(MM[,v2]))
    a1 <- a2 <- matrix(0, nrow = 3, ncol = ncol(MM))
    a1[, ind1[1]] <- 1
    a1[, ind1[2]] <- s2
    a2[, ind2[1]] <- 1
    a2[, ind2[2]] <- s1
    eff1 <- a1 %*% obj$coef
    if(is.null(varcov)){
        varcov <- vcov(obj)
    }
    se.eff1 <- sqrt(diag(a1 %*% varcov %*% t(a1)))
    low1 <- eff1 - qt(0.975, obj$df.residual) * se.eff1
    up1 <- eff1 + qt(0.975, obj$df.residual) * se.eff1
    eff2 <- a2 %*% obj$coef
    se.eff2 <- sqrt(diag(a2 %*% varcov %*% t(a2)))
    low2 <- eff2 - qt(0.975, obj$df.residual) * se.eff2
    up2 <- eff2 + qt(0.975, obj$df.residual) * se.eff2
    if (!plot.type %in% c("pdf", "png", "eps", "screen")) {
        print("plot type must be one of - pdf, png or eps")
    }
    else {
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v1, ".pdf", sep = ""),
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v1, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            old.psopts <- ps.options()
            setEPS()
            postscript(paste(name.stem, "_", v1, ".eps", sep = ""))
        }
        if (plot.type == "screen") {
            oldpar <- par()
            par(mfrow = c(1, 2))
        }
        plot(s2, eff1, type = "n", ylim = range(c(low1, up1)), axes=FALSE,
            xlab = ifelse(is.null(xlab), toupper(v2), xlab[1]), ylab = ifelse(is.null(ylab), paste("Conditional Effect of ",
                toupper(v1), " | ", toupper(v2), sep = ""), ylab[1]))
            axis(1, at=s2, labels=c("Mean-SD", "Mean", "Mean+SD"))
            axis(2)
            box()

        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        points(s2, eff1, pch=16)
        segments(s2, low1, s2, up1)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "pdf") {
            pdf(paste(name.stem, "_", v2, ".pdf", sep = ""),
                height = 6, width = 6)
        }
        if (plot.type == "png") {
            png(paste(name.stem, "_", v2, ".png", sep = ""))
        }
        if (plot.type == "eps") {
            postscript(paste(name.stem, "_", v2, ".eps", sep = ""))
        }
        plot(s1, eff2, type = "n", ylim = range(c(low2, up2)),
            axes=FALSE,
            xlab = ifelse(is.null(xlab), toupper(v1), xlab[2]),
			ylab = ifelse(is.null(ylab), paste("Conditional Effect of ",
                toupper(v2), " | ", toupper(v1), sep = ""), ylab[2]))
            axis(1, at=s1, labels=c("Mean-SD", "Mean", "Mean+SD"))
            axis(2)
            box()
        if (par()$usr[3] < 0 & par()$usr[4] > 0) {
            abline(h = 0, col = "gray50")
        }
        points(s1, eff2, pch=16)
        segments(s1, low2, s1, up2)
        if (plot.type != "screen") {
            dev.off()
        }
        if (plot.type == "eps") {
            ps.options <- old.psopts
        }
        if (plot.type == "screen") {
            par <- oldpar
        }
    }
}

##' Print method for glmChange objects
##' 
##' @description Print method for object of class \code{glmc2}. 
##' 
##' @param x Object of class \code{glmc2}
##' @param ... Currently unimplemented. 
##' @export
##' @method print glmc2
print.glmc2 <- function(x, ...){
    print.default(x$res)
    invisible(x)
}

##' Print method for intQualQuant objects. 
##' 
##' @description Print method for objects of class \code{iqq} calculated
##' with the \code{intQualQuant} function.
##' 
##' @param x Object of class \code{iqq}
##' @param ... Currently unimplmeneted
##' @export 
##' @method print iqq
print.iqq <- function(x, ...){
    cat("Conditional Effect of ", x$mainvar, " given ", x$givenvar, "\n")
    printCoefmat(x$out, digits=4)
    cat("\n")
    invisible(x)
}

central <- function(x, ...){
    UseMethod("central")
}

##' @method central factor
central.factor  <- function(x, ...){
    tab <- table(droplevels(x))
    res <- x[1]
    res[1] <- names(tab)[which.max(tab)]
    return(res)
}

##' @method central numeric
central.numeric <- function(x, type=c("median", "mean"), ...){
    type <- match.arg(type)
    if(type == "median"){
        res <- median(x, na.rm=TRUE, ...)
    }
    else{
        res <- mean(x, na.rm=TRUE, ...)
    }
    return(res)
}


#' Calculate Cross-Derivative and its Variability
#' 
#' Calculates the cross-derivative required to evaluate interactions in
#' logistic/probit regression models.
#' 
#' The function calculates the second difference as (Pr(Y=1|x1=max, x2=max) -
#' Pr(Y=1|x1=min, x2=max)) - (Pr(Y=1|x1=max, x2=min) - Pr(Y=1|x1=min, x2=min)).
#' The function uses a parametric bootstrap to calculate the sampling
#' distribution of the second difference.
#' 
#' @aliases secondDiff 
#' 
#' @param obj An object of class \code{glm} that will be used to find the
#' cross-derivative.
#' @param vars A vector of two variables to be used in calculating the
#' derivative.
#' @param data A data frame.
#' @param method Indicate whether you want to use average marginal effects
#' (AME) or marginal effects at representative values (MER).
#' @param vals A named list of length 2 where each element gives the minimum
#' and maximum values used in the calculation.
#' @param typical A named vector of values at which to hold variables constant.
#' @return A list with two elements:
#' 
#' \item{ave}{The average second difference in each iteration of the
#' bootstrap.} \item{ind}{If \code{type == 'AME'}, \code{ind} is returned with
#' the second difference and measures of uncertainty for each individual
#' observation in the original dataset} \item{probs}{If \code{type == 'MER'},
#' \code{probs} is returned with the full matrix of simulated predicted
#' probabilities for the four conditions.}
#'
#' @export
#' @author Dave Armstrong
secondDiff <- function(obj, vars, data, method=c("AME", "MER"), vals = NULL, typical = NULL){
    disc <- sapply(vars, function(x)is.factor(data[[x]]))
    meth <- match.arg(method)
    mf <- model.frame(obj)
    b <- coef(obj)
    if(is.null(vals)){
      if(!disc[1]){
          min1 <- min(data[[vars[1]]], na.rm=TRUE)
          max1 <- max(data[[vars[1]]], na.rm=TRUE)
      }
      if(disc[1]){
          cn1 <- paste(vars[1], levels(droplevels(mf[, vars[1]])), sep="")
          b1 <- b[cn1]
          b1 <- ifelse(is.na(b1), 0, b1)
          names(b1) <- levels(droplevels(mf[, vars[1]]))
          min1 <- max1 <- mf[1,vars[1]]
          min1[1] <- names(b1)[which.min(b1)]
          max1[1] <- names(b1)[which.max(b1)]
      }
      if(!disc[2]){
          min2 <- min(data[[vars[2]]], na.rm=TRUE)
          max2 <- max(data[[vars[2]]], na.rm=TRUE)
      }
      if(disc[2]){
          cn2 <- paste(vars[2], levels(droplevels(mf[, vars[2]])), sep="")
          b2 <- b[cn2]
          b2 <- ifelse(is.na(b2), 0, b2)
          names(b2) <- levels(droplevels(mf[, vars[2]]))
          min2 <- max2 <- mf[1,vars[2]]
          min2[1] <- names(b2)[which.min(b2)]
          max2[1] <- names(b2)[which.max(b2)]
      }
    }else{
      if(disc[1]){
        l1 <- levels(mf[, vars[1]])
        w1 <- which(l1 == vals[[vars[1]]][1])
        w2 <- which(l1 == vals[[vars[1]]][2])
        if(length(w1) == 0){
          stop(glue("{vals[[vars[1]]][1]} not a level of {vars[1]}\n"))
        }
        if(length(w2) == 0){
          stop(glue("{vals[[vars[1]]][2]} not a level of {vars[1]}\n"))
        }
        v1 <- factor(w1, levels=1:length(l1), labels=l1)
        v2 <- factor(w2, levels=1:length(l1), labels=l1)
        min1 <- v1
        max2 <- v2
        
      }else{
        min1 <- vals[[vars[1]]][1]
        max1 <- vals[[vars[1]]][2]
      }
      if(disc[2]){
        l1 <- levels(mf[, vars[2]])
        w1 <- which(l1 == vals[[vars[2]]][1])
        w2 <- which(l1 == vals[[vars[2]]][2])
        if(length(w1) == 0){
          stop(glue("{vals[[vars[2]]][1]} not a level of {vars[2]}\n"))
        }
        if(length(w2) == 0){
          stop(glue("{vals[[vars2]][2]} not a level of {vars[2]}\n"))
        }
        v1 <- factor(w1, levels=1:length(l1), labels=l1)
        v2 <- factor(w2, levels=1:length(l1), labels=l1)
        min2 <- v1
        max2 <- v2
        
      }else{
        min2 <- vals[[vars[2]]][1]
        max2 <- vals[[vars[2]]][2]
      }
    }
    b <- mvrnorm(1500, coef(obj), vcov(obj))

    if(meth == "AME"){
        dat1 <- dat2 <- dat3 <- dat4 <- data
        # dat1 = (L, L)
        # dat2 = (H, L)
        # dat3 = (L, H)
        # dat4 = (H, H)
        dat1[[vars[1]]] <- min1
        dat1[[vars[2]]] <- min2
        dat2[[vars[1]]] <- max1
        dat2[[vars[2]]] <- min2
        dat3[[vars[1]]] <- min1
        dat3[[vars[2]]] <- max2
        dat4[[vars[1]]] <- max1
        dat4[[vars[2]]] <- max2
    modmats <- list()
    modmats[[1]] <- model.matrix(formula(obj), dat1)
    modmats[[2]] <- model.matrix(formula(obj), dat2)
    modmats[[3]] <- model.matrix(formula(obj), dat3)
    modmats[[4]] <- model.matrix(formula(obj), dat4)
    probs <- lapply(modmats, function(x)family(obj)$linkinv(x%*%t(b)))
    D <- c(1,-1,-1,1)
    secdiff <- sapply(1:1500, function(x)(cbind(probs[[1]][,x], probs[[2]][,x], probs[[3]][,x], probs[[4]][,x]) %*%D))
    avesecdiff <- colMeans(secdiff)
    indsecdiff <- rowMeans(secdiff)
    indp <- apply(secdiff, 1, function(x)mean(x > 0))
    indp <- ifelse(indp > .5, 1-indp, indp)
    indci <- t(apply(secdiff, 1, quantile, c(.025,.975)))
    ret <- list(ave = avesecdiff, ind=data.frame(secddiff = indsecdiff, pval = indp, lower=indci[,1], upper=indci[,2], stringsAsFactors=TRUE))
  }
  else{
    v <- all.vars(formula(obj))[-1]
    rn <- v
    tmp.df <- do.call(data.frame, c(lapply(v, function(x)central(data[[x]])), stringsAsFactors=TRUE))
    names(tmp.df) <- v
    if(!is.null(typical)){
        for(i in 1:length(typical)){
            tmp.df[[names(typical)[[i]]]] <- typical[[i]]
        }
    }
    tmp.df <- tmp.df[c(1,1,1,1), ]
    tmp.df[[vars[1]]][1] <- min1
    tmp.df[[vars[1]]][2] <- max1
    tmp.df[[vars[1]]][3] <- min1
    tmp.df[[vars[1]]][4] <- max1
    tmp.df[[vars[2]]][1] <- min2
    tmp.df[[vars[2]]][2] <- min2
    tmp.df[[vars[2]]][3] <- max2
    tmp.df[[vars[2]]][4] <- max2
    f <- formula(obj)
    f[[2]] <- NULL
    modmat <- model.matrix(f, data=tmp.df)
    preds <- t(family(obj)$linkinv(modmat%*%t(b)))
    D <- c(1,-1,-1,1)
    secdiff <- (preds %*% D)
    ret <- list(ave = secdiff, probs=preds)
  }
    class(ret) <- "secdiff"
    return(ret)
}

##' Summary for Second Difference Objects
##' 
##' @description Summary method for objects of class \code{secdiff}. 
##' 
##' @param object An object of class \code{secdiff}
##' @param level Confidence level for the confidence intervals
##' @param digits Number of digits to print
##' @param ... Other arguments to be passed down to \code{summary}. 
##' 
##' @method summary secdiff
##' @export
summary.secdiff <- function(object, ..., level=0.95, digits=3){
  ll <- (1-level)/2
  ul <- 1-ll
  type <- ifelse("ind" %in% names(object), "Average Marginal Effect", "Marginal Effect at Typical Values")
  cat("Second Difference Using the", type, "Approach\n\n")
  s <- c(mean(object$ave), quantile(object$ave, c(ll, ul)))
  s <- sprintf(glue("%.{digits}f"), s)
  if(type == "Average Marginal Effect"){
    cat("Overall: \n")
  }
  cat("Average Second Difference: ", s[1], ", ", level*100, "% CI: (", s[2], ",", s[3], ")\n\n", sep="")
  if(type == "Average Marginal Effect"){
    cat("Individual:\n")
  sigpos <- sum(object$ind$lower > 0)
  signeg <- sum(object$ind$upper < 0)
  insig <- sum(object$ind$lower <=0 & object$ind$upper >=0)
  cat("Significant Negative Individual Second Differences:", signeg, "\n")
  cat("Significant Positive Individual Second Differences:", sigpos, "\n")
  cat("Inignificant Individual Second Differences:", insig, "\n")
  }
}


##' Plotting Method for Second Difference Objects
##' 
##' @description Plots the results of the \code{secondDiff} function.
##' 
##' @param x An object of class \code{secdiff} 
##' @param level The confidence level of the confidence interval(s)
##' @param ... Other arguments to be passed down, currently not implemented. 
##' 
##' @method plot secdiff
##' @export
plot.secdiff <- function(x, level=.95, ...){
  ll <- (1-level)/2
  ul <- 1-ll
  type <- ifelse("ind" %in% names(x), "Average Marginal Effect", "Marginal Effect at Typical Values")
  s <- c(mean(x$ave), quantile(x$ave, c(ll, ul)))
  ave.df <- data.frame(difference = s[1], lower=s[2], upper=s[3], stringsAsFactors=TRUE)
  if(!("ind" %in% names(x))){
    ave.df$x <- factor(1, levels=1, labels="Average")
    g <- ggplot(ave.df) + 
      geom_point(aes_string(y="difference", x="x")) + 
      geom_segment(aes_string(x="x", xend="x", y="lower", yend="upper")) + 
      theme_bw() + 
      labs(x="", y="Second Difference") + 
      ggtitle(glue("Plot of Second Differences\n{type} Approach"))
  }
    if("ind" %in% names(x)){
    ind.df <- x$ind
    ind.df <- ind.df[order(ind.df$secddiff), ] 
    ind.df$obs <- 1:nrow(ind.df)
    d <- abs(ave.df$difference - ind.df$secddiff)
    w <- which.min(d)
    ave.df <- data.frame(secddiff = s[1], lower=s[2], upper=s[3], obs=w, stringsAsFactors=TRUE)
    g <- ggplot(ind.df) + 
      geom_point(aes_string(y="secddiff", x="obs"), size=.5) + 
      geom_segment(aes_string(x="obs", xend="obs", y="lower", yend="upper"), 
                   size=.25, col="gray75", alpha=.5) + 
      geom_hline(aes(yintercept=0), lty=2, size=.75) + 
      geom_point(data=ave.df, aes_string(y="secddiff", x="obs"), size=.75, col="red") +
      geom_segment(data=ave.df, aes_string(x="obs", xend="obs", y="lower", yend="upper"), 
                   size=1, col="red", alpha=.5) +
      theme_bw() + 
      labs(x="", y="Second Difference") + 
      ggtitle(glue("Plot of Individual Second Differences\n{type} Approach"))
      }
}




#' Plot Effects of Removing Outliers
#' 
#' Plots the effect of a variable sequentially removing outlying observations.
#' 
#' 
#' @param obj An object of class \code{glm} that will be used to plot the
#' outlier removed lines.
#' @param var A character string giving the name of the variable to be used.
#' @param data A data frame.
#' @param stat Which statistic to use to evaluate \sQuote{outlyingness}.
#' @param nOut Number of outliers to be removed.
#' @param whichOut If not \code{NULL}, a vector of observation numbers to be
#' removed manually, rather than using \code{stat}.
#' @param cumulative Logical indicating whether the outliers should be removed
#' cumulatively, or each one in turn, replacing the other outlying
#' observations.
#' 
#' @export
#' 
#' @author Dave Armstrong
outEff <- function(obj, var, data, stat =c("cooksD", "hat", "deviance", "pearson"), nOut = 10, whichOut=NULL, cumulative=FALSE){
    stat <- match.arg(stat)
    if(is.null(whichOut)){
        vec <- switch(stat,
            cooksD = cooks.distance(obj),
            hat = hatvalues(obj),
            deviance = residuals(obj, type="deviance"),
            pearson = residuals(obj, type="pearson"))
        tmp <- abs(vec)
        rnk <- rank(-tmp)
        w <- which(rnk %in% 1:nOut)
    }
    else{
        w <- whichOut
        nOut <- length(w)
    }
    mf <- model.frame(obj)
    eff.list <- fits <- list()
    eff.list[[1]] <- effect(var, obj, xlevels=25)
    sigs <- rep(0, nOut)
    k <- 2
    for(i in 1:nOut){
        if(cumulative){
            outs <- w[1:i]
        }
        else{
            outs <- w[i]
        }
        tmp <- data[-outs, ]
        tmpObj <- update(obj, data=tmp)
        s <- Anova(tmpObj)
        pval <- s[which(rownames(s) == var), 3]
        sigs[i] <- ifelse(pval < 0.05, 1, 0)
        eff.list[[k]] <- effect(var, tmpObj, xlevels=25)
        fits[[k]] <- eff.list[[k]]$transformation$inverse(eff.list[[k]]$fit)
        k <- k+1
    }
    fits <- do.call(cbind, fits)
    yrg <- c(min(c(fits)), max(c(fits)))
    if(is.numeric(mf[[var]]) & length(unique(mf[[var]])) > 2){
        plot(eff.list[[1]]$x[[var]], fits[,1], type="n", ylim=yrg,
            xlab = var, ylab = "Fitted Values")
        cols <- c("gray65", "red")
        for(i in 2:ncol(fits)){
            lines(eff.list[[1]]$x[[var]], fits[,i], col=cols[(sigs[i]+1)])
        }
        lines(eff.list[[1]]$x[[var]], fits[,1], col="black", lwd=1.5)
    }
    else{
        ins <- strwidth(eff.list[[1]]$x[[var]], units="inches")
        ins <- max(ins)
        pmai <- par()$mai
        pmai[1] <- .25+ins
        par(mai = pmai)
        plot(1:length(eff.list[[1]]$x[[var]]), fits[,1], type="n", ylim=yrg,
            xlab = "", ylab = "Fitted Values", axes=FALSE)
        axis(2)
        axis(1, at=1:length(eff.list[[1]]$x[[var]]), labels=eff.list[[1]]$x[[var]], las=2)
        box()
                cols <- c("gray65", "red")
                for(i in 2:ncol(fits)){
                    points(jitter(1:length(eff.list[[1]]$x[[var]]), 1), fits[,i], col=cols[(sigs[i]+1)])
                }
                points(1:length(eff.list[[1]]$x[[var]]), fits[,1], col="black", pch=16)
            }
    }



#' Plot First Differences from Ordinal DV Model
#' 
#' Takes the output from \code{\link{ordChange}} and turns it into a plot.
#' 
#' 
#' @param ordc The output from \code{ordChange}.
#' @param plot Logical indicating whether a plot (if \code{TRUE}) or data (if
#' \code{FALSE}) should be returned.
#' @return Either a \code{lattice} plot or a \code{data.frame} depending on the
#' specification of the \code{plot} argument.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @examples
#' library(MASS)
#' data(france)
#' polr.mod <- polr(vote ~ age + male + retnat + lrself, data=france)
#' typical.france <- data.frame(
#' 	age = 35, 
#' 	retnat = factor(1, levels=1:3, labels=levels(france$retnat)), 
#' 	stringsAsFactors=TRUE)
#' oc.res <- ordChange(polr.mod, data=france, typical.dat=typical.france, sim=TRUE)	
#' oc2plot(oc.res)
oc2plot <- function(ordc, plot=TRUE){
    tmpdat <- data.frame(
        var = rep(rownames(ordc$diffs$mean), ncol(ordc$diffs$mean)),
        lev = rep(colnames(ordc$diffs$mean), each = nrow(ordc$diffs$mean)),
        mean = c(ordc$diffs$mean),
        lower=c(ordc$diffs$lower),
        upper = c(ordc$diffs$upper), 
        stringsAsFactors=TRUE)
    p1 <- xyplot(mean ~ lev | var, data=tmpdat,
        xlab = "", ylab = "Predicted Change in Pr(y=m)",
        lower = tmpdat$lower, upper=tmpdat$upper,
        panel = function(x,y,subscripts, lower, upper,...){
            panel.abline(h=0, lty=2)
            panel.arrows(x, lower[subscripts], x, upper[subscripts], angle=90, length=.05, code=3)
            panel.points(x,y,pch=16, cex=.75, col="black")
        }, prepanel=prepanel.ci)
    if(plot){
        return(p1)
    }
    else{
        return(tmpdat)
    }
}



#' Plog Probabilities by Group
#' 
#' Plots predicted probabilities by value of the dependent variable for
#' proportional odds logistic regression and multinomial logistic regression
#' models
#' 
#' Plots the predicted probabilities by value of the dependent variable.  Each
#' panel only includes the observations that had the identified value of the
#' dependent variable.
#' 
#' @aliases probgroup probgroup.polr probgroup.multinom
#' @param obj Object of class \code{polr} or \code{multinom} where appropraite.
#' @param ... Currently not implemented.
#' @return A plot.
#' 
#' @export
#' 
#' @author Dave Armstrong
probgroup <- function(obj, ...){
    UseMethod("probgroup")
}

##' @method probgroup polr
##' @export
probgroup.polr <- function(obj, ...){
    pr <- predict(obj, type="probs")
    y <- model.response(model.frame(obj))
    pr2 <- pr[cbind(1:nrow(pr), as.numeric(y))]
    tmpdat <- data.frame(probs=c(pr2),
        y = factor(as.numeric(y), labels=obj$lev), 
        stringsAsFactors=TRUE)
    ly <- length(table(y))
    return(histogram(~probs | y, data=tmpdat, col="gray75", ylim = c(0, 100),
        layout=c(ly,1), xlab="Predicted Probabilities",
        par.settings=list(strip.background=list(col="white"))))
}

##' @method probgroup multinom
##' @export
probgroup.multinom <- function(obj, ...){
    pr <- predict(obj, type="probs")
    y <- model.response(model.frame(obj))
    pr2 <- pr[cbind(1:nrow(pr), as.numeric(y))]

    tmpdat <- data.frame(probs=c(pr2),
        y = y, stringsAsFactors=TRUE)
    ly <- length(table(y))
    return(histogram(~probs | y, data=tmpdat, col="gray75", ylim = c(0, 100),
        layout=c(ly,1), xlab="Predicted Probabilities",
        par.settings=list(strip.background=list(col="white"))))
}



#' Scalar Measures of Fit for Poisson GLMs Models
#' 
#' Calculates scalar measures of fit for models with count dependent variables
#' along the lines described in Long (1997) and Long and Freese (2005).
#' 
#' \code{poisfit} calculates scalar measures of fit (many of which are
#' pseudo-R-squared measures) to describe how well a model fits data with a
#' count dependent variable.
#' 
#' @param obj A model of class \code{glm} with \code{family=poisson}.
#' @return A named vector of scalar measures of fit
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @references Long, J.S.  1997.  Regression Models for Categorical and Limited
#' Dependent Variables.  Thousand Oaks, CA: Sage.
#' 
#' Long, J.S. and J. Freese.  2005.  Regression Models for Categorical Outcomes
#' Using Stata, 2nd ed.  College Station, TX: Stata Press.
poisfit <- function(obj){
	y <- model.response(model.frame(obj))
	if(family(obj)$family == "poisson"){
	  nullMod <- glm(y ~ 1, family=poisson)
	}else if(inherits(obj, "negbin")){
	  nullMod <- glm.nb(y ~ 1)
	} else{
	  stop("Model must be either poisson or negative binomial")
	}
	llNull <- logLik(nullMod)
	r2_mcf <- 1-(logLik(obj)/llNull)
	r2_mcfa <- 1-((logLik(obj) - obj$rank)/llNull)
	g2 <- -2*(llNull - logLik(obj))
	r2_ml <- 1-exp(-g2/length(y))
    r2cu <- (r2_ml)/(1-exp(2*llNull/length(y)))
	res <- matrix(nrow=6, ncol=2)
	colnames(res) <- c("Estimate", "p-value")
	rownames(res) <- c("GOF (Pearson)", "GOF (Deviance)", "ML R2", "McFadden R2", "McFadden R2 (Adj)", "Cragg-Uhler(Nagelkerke) R2")
	gof <- poisGOF(obj)
	res[1:2,1] <- as.numeric(as.character(gof[,1]))
	res[1:2,2] <- as.numeric(as.character(gof[,2]))
	res[3,1] <- r2_ml
	res[4,1] <- r2_mcf
	res[5,1] <- r2_mcfa
	res[6,1] <- r2cu
	if(inherits(obj, "negbin")){res <- res[-(1:2), ]}
	pres <- matrix(apply(res, 2, function(x)sprintf("%.3f", x)), ncol=2)
	dimnames(pres) <- dimnames(res)
	print(pres, quote=FALSE)
	invisible(res)
}



#' Make Hypothetical Predictions for Survey Data
#' 
#' Calculates survival probabilities for hypothetical data.
#' 
#' 
#' @param l A named list where variable names are the names and values of the
#' variables are the values.  Combinations will be made with
#' \code{expand.grid}.
#' @param obj A model object estimated with \code{survreg}.
#' @param ... currently not implemented.
#' @return A data frame.
#' 
#' @export
#' 
#' @author Dave Armstrong
makeHypSurv <- function(l,obj, ...){
    tmp <- do.call(expand.grid, l)
    p <- seq(.99,0,by=-.01)
    preds <- predict(obj, newdata=tmp,
    type="quantile", p=p)
    plot.data <- data.frame(
        p.fail = rep(p, each=nrow(preds)),
        time = c(preds), stringsAsFactors=TRUE)
    plot.data <- cbind(plot.data, tmp[rep(1:3, nrow(plot.data)/3), ])
    rownames(plot.data) <- NULL
    return(plot.data)
}



#' Lag a time-series cross-sectional variables
#' 
#' Lags (or leads) a variable in a time-series corss-sectional dataset.
#' 
#' 
#' @param dat A data frame.
#' @param x A string identifying variable to be lagged.
#' @param id A string identifying the name of the cross-sectional identifier.
#' @param time A string identifying the name of the time variable.
#' @param lagLength The length of the lag, use negative values for leading
#' variables.
#' @return A vector giving the lagged values of \code{x}.
#' 
#' @export
#' 
#' @author Dave Armstrong
tscslag <- function(dat, x, id, time, lagLength=1){
	obs <- apply(dat[, c(id, time)], 1,
	    paste, collapse=".")
	tm1 <- dat[[time]] - lagLength
	lagobs <- apply(cbind(dat[[id]], tm1),
	    1, paste, collapse=".")
	lagx <- dat[match(lagobs, obs), x]
	lagx
}



#' Test Transformations and Polynomials in Non-linear Models
#' 
#' Tests for model improvements for non-linear transformations and polynomials
#' with Clarke's (2007) distribution-free test for non-nested models.
#' 
#' Three hypotheses are tested with this function.  The first is whether the
#' original specification is preferred to the power transformation.  The second
#' is whether the original specification is preferred to the polynomial model.
#' The third is whether the power transformation is preferred to the polynomial
#' model.  All tests are done with the Clarke test.
#' 
#' @aliases testNL testNL.glm testNL.lm
#' @param obj Object of a supported class in which non-linear functional forms
#' will be tested.
#' @param var String giving name of variable to be tested.
#' @param transPower The power used in the transformation.  For transformations
#' in the range (-0.01, 0.01), the log transformation is used.
#' @param polyOrder The order of the polynomial to be used.
#' @param plot Logical indicating whether the effects should be plotted
#' @param ... Currently not implemented.
#' @return A plot or a data frame giving the results of the tests identified
#' above.
#' @author Dave Armstrong
#' 
#' @export
#' 
#' @references Kevin Clarke.  2007.  "A Simple Distribution-Free Test for
#' Nonnested Hypotheses."  \emph{Political Analysis} 15(3): 347--363.
testNL <- function(obj, var, transPower, polyOrder, plot=FALSE, ...){
  UseMethod("testNL")
}

##' Power Transformation Function
##' 
##' @description Power transformation function that treats everything with
##' absolute power transform < .01 as the log transform.  
##' @param x Vector of values to be transformed
##' @param transPower The power of the transformation
##' @return A vector of transformed values
##' @export
powerTrans <- function(x, transPower){
  if(abs(transPower) > .01){
    x <- I(x^transPower)
  }  else {
    x <- log(x)
  }
  x
}

##' @rdname testNL
##' @export
##' @method testNL glm
testNL.glm <- function(obj, var, transPower, polyOrder, plot=FALSE, ...){
  res <- data.frame(Model1 = c("Original", "Original", "Power Transform"),
                    Model2 = c("Power Transform", "Polynomial", "Polynomial"),
                    pval = NA,
                    preferred = NA, stringsAsFactors=TRUE)

  form1 <- glue("~ . -{var} + powerTrans({var}, {transPower})")
  form2 <- glue("~ . -{var} + poly({var}, {polyOrder}, raw=TRUE)")
  o1 <- update(obj, form1)
  t1 <- clarke_test(obj, o1)
  b <- min(t1$stat, t1$nobs - t1$stat)
  p <- 2 * pbinom(b, t1$nobs, 0.5)
  pref <- ifelse(t1$stat > t1$nobs - t1$stat, 1, 2)
  res$pval[1] <- p
  if(p < .05 & pref == 1){
    res$preferred[1] <- "Original"
  }
  if(p < .05 & pref == 2){
    res$preferred[1] <- "Power Transform"
  }
  if(p >= .05){
    res$preferred[1] <- "Neither"
  }

  o2 <- update(obj, form2)
  t2 <- clarke_test(obj, o2)
  b <- min(t2$stat, t2$nobs - t2$stat)
  p <- 2 * pbinom(b, t2$nobs, 0.5)
  pref <- ifelse(t2$stat > t2$nobs - t2$stat, 1, 2)
  res$pval[2] <- p
  if(p < .05 & pref == 1){
    res$preferred[2] <- "Original"
  }
  if(p < .05 & pref == 2){
    res$preferred[2] <- "Polynomial"
  }
  if(p >= .05){
    res$preferred[2] <- "Neither"
  }

  t3 <- clarke_test(o1, o2)
  b <- min(t3$stat, t3$nobs - t3$stat)
  p <- 2 * pbinom(b, t3$nobs, 0.5)
  pref <- ifelse(t3$stat > t3$nobs - t3$stat, 1, 2)
  res$pval[3] <- p
  if(p < .05 & pref == 1){
    res$preferred[3] <- "Power Transform"
  }
  if(p < .05 & pref == 2){
    res$preferred[3] <- "Polynomial"
  }
  if(p >= .05){
    res$preferred[3] <- "Neither"
  }
  res$pval <- sprintf("%.3f", res$pval)
  xl <- list(25)
  names(xl)[1] <- var
  if(!plot){
    return(res)
  }
  if(plot){
    e1 <- Effect(var, obj, xlevels=xl)
    se1 <- do.call(data.frame, c(summary(e1)[c("effect", "lower", "upper")], stringsAsFactors=TRUE))
    se1$model <- factor(1, levels=1:3,
                        labels=c("Original", "Power", "Polynomial"))
    se1$x <- e1$x[[var]]

    e2 <- Effect(var, o1, xlevels=xl)
    se2 <- do.call(data.frame, c(summary(e2)[c("effect", "lower", "upper")], stringsAsFactors=TRUE))
    se2$model <- factor(2, levels=1:3,
                        labels=c("Original", "Power", "Polynomial"))
    se2$x <- e2$x[[var]]

    e3 <- Effect(var, o2, xlevels=xl)
    se3 <- do.call(data.frame, c(summary(e3)[c("effect", "lower", "upper")], stringsAsFactors=TRUE))
    se3$model <- factor(3, levels=1:3,
                        labels=c("Original", "Power", "Polynomial"))
    se3$x <- e3$x[[var]]
    plotData <- rbind(se1, rbind(se2, se3))
    ggplot(plotData, aes_string(x="x", y="effect", ymin = "lower", 
            ymax="upper", fill="model")) + 
      geom_ribbon(colour=NA, alpha=.25) + 
      geom_line(aes_string(colour="model")) +
      theme_bw() +  
      theme(legend.position="bottom") + 
      labs(x=var)

  }
}

##' @rdname testNL
##' @export
##' @method testNL lm
testNL.lm <- testNL.glm



#' Plot Effects from Firth Logit
#' 
#' Plots the effect of a variable in a model estimated with Firth Logit.
#' 
#' The \code{effect.logistf} function calculates the effect (predicted
#' probabilities) of a variable in a Firth logit model estimated with the
#' \code{logistf} function.  The function estimates the analogous glm.
#' It then replaces the coefficient vector in that model object with the Frith
#' logit coefficients.  It also puts the variance-covariance matrix from the
#' Firth logit in the model object and uses a custom extractor function in the
#' \code{Effect} function to extract that variance-covariance matrix rather
#' than the one usually extracted with \code{vcov}.  Note that variability and
#' confidence intervals for the effects will not be calculated using profile
#' likelihood as they are in the Firth logit, but will be calculated using the
#' appropriate variance-covariance matrix.
#' 
#' @param var A character string giving the name of the variable whose effect
#' is to be generated.
#' @param obj An object of class \code{logistf}.
#' @param data A data frame.
#' @param ... Other arguments to be passed down to the \code{\link{Effect}}
#' function.
#' @return An object of class \code{eff} that can be used with other functions
#' from the \code{effects} package.
#' 
#' @export
#' 
#' @author Dave Armstrong
effect_logistf <- function(var, obj, data, ...){
  v <- function(obj, complete=FALSE)return(obj$var)
  form <- as.character(obj$formula)
  form <- as.formula(paste(form[2], form[3], sep=form[1]))
  tmp.mod <- glm(form, data=data, family=binomial)
  tmp.mod$coefficients <- as.vector(obj$coefficients)
  tmp.mod$var <- obj$var
  e <- Effect(var, tmp.mod, vcov.=v, ...)
  return(e)
}

