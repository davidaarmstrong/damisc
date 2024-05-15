#### Utility Functions defined by Fox in the 'car' package (as of 6-25-2012)
#### These are used in the crTest function, 
# citation: 
# John Fox and Sanford Weisberg (2011). An {R} Companion to
# Applied Regression, Second Edition. Thousand Oaks CA: Sage.
# URL: http://socserv.socsci.mcmaster.ca/jfox/Books/Companion

#' Identify if model has intercept
#' 
#' @param model A model object
#' @param ... Other arguments passed down, currently not implemented
#' @export
has.intercept <- function (model, ...) {
	UseMethod("has.intercept")
}

#' @method has.intercept default
#' @export
has.intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

#' Retrieve term names from model
#' 
#' @param model A model object
#' @param ... Other arguments passed down, currently not implemented
#' @export
term.names <- function (model, ...) {
	UseMethod("term.names")
}

#' @method term.names default
#' @export
term.names.default <- function (model, ...) {
	term.names <- labels(terms(model))
	if (has.intercept(model)) c("(Intercept)", term.names)
	else term.names
}

#' Retrieve predictor names from model
#' 
#' @param model A model object
#' @param ... Other arguments passed down, currently not implemented
#' @export
predictor.names <- function(model, ...) {
	UseMethod("predictor.names")
}

#' @method predictor.names default
#' @export
predictor.names.default <- function(model, ...){
	predictors <- attr(terms(model), "variables")
	as.character(predictors[3:length(predictors)])
}

#' Retrieve response name from model
#' 
#' @param model A model object
#' @param ... Other arguments passed down, currently not implemented
#' @export
responseName <- function (model, ...) {
	UseMethod("responseName")
}

#' @method responseName default
#' @export
responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

#' Retrieve model response
#' 
#' @param model A model object
#' @param ... Other arguments passed down, currently not implemented
#' @export
response <- function(model, ...) {
	UseMethod("response")
}

#' @method response default
#' @export
response.default <- function (model, ...) model.response(model.frame(model))

#' Is coefficient aliased
#' 
#' Identifies whether there are aliased coefficients in model
#' @param model A model object
#' @importFrom stats alias
#' @export
is.aliased <- function(model){
	!is.null(alias(model)$Complete)
}

#' Retrieve term degrees of freedom from model
#' 
#' @param model A model object
#' @param term A model term
#' @param ... Other arguments passed down, currently not implemented
#' @export
df.terms <- function(model, term, ...){
	UseMethod("df.terms")
}

#' @method df.terms default
#' @export 
df.terms.default <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}
