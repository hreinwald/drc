#' @title Confidence Intervals for model parameters
#'
#' @description
#' Computes confidence intervals for one or more parameters in a model of class
#' 'drc'.
#'
#' @param object a model object of class 'drc'.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level the confidence level required.
#' @param pool logical. If TRUE curves are pooled. Otherwise they are not. This
#'   argument only works for models with independently fitted curves as
#'   specified in \code{\link{drm}}.
#' @param ... additional argument(s) for methods. Not used.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence
#'   limits for each parameter. These will be labelled as (1-level)/2 and
#'   1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' ## Confidence intervals for all parameters
#' confint(ryegrass.m1)
#'
#' ## Confidence interval for a single parameter
#' confint(ryegrass.m1, "e")
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"confint.drc" <- function(object, parm, level = 0.95, pool = TRUE, ...)
#"confint.drc" <- function(object, parm, level = 0.95, type = "t", pool = TRUE, ...)
{
    ## Matching parameter names
    if (!missing(parm))
    {
        matchVec <- object$"parNames"[[2]] %in% parm
        if (!any(matchVec)) {stop("The 'parm' argument does not match an actual parameter")}
    } else {
        matchVec <- rep(TRUE, length(object$"parNames"[[2]]))
    }

    ## Constructing matrix of confidence intervals    
    confint.basic(summary(object, pool = pool)$"coefficients"[matchVec, 1:2, drop = FALSE], 
    level, object$"type", df.residual(object))  
    
#    ## Retrieving estimates and estimated standard errors
#    estMat <- summary(object, pool = pool)$"coefficients"[matchVec, 1:2, drop = FALSE]
#
#    ## Constructing matrix of confidence intervals
#    confMat <- matrix(0, dim(estMat)[1], 2)
#
#    alphah <- (1 - level)/2 
#    if (type == "u") {two <- qnorm(1 - alphah)}
#    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}
#    confMat[, 1] <- estMat[, 1] -  two * estMat[, 2]
#    confMat[, 2] <- estMat[, 1] +  two * estMat[, 2]
#
#    ## Formatting matrix
#    colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " ") )
#    rownames(confMat) <- rownames(estMat)
#
#    return(confMat)  
}

## Defining basic function for providing confidence intervals
"confint.basic" <- function(estMat, level, intType, dfres, formatting = TRUE)
{
    alphah <- (1 - level)/2 
#    if (type == "u") {two <- qnorm(1 - alphah)}
#    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}    
    tailPercentile <- switch(intType, 
    binomial = qnorm(1 - alphah), 
    continuous = qt(1 - alphah, dfres),
    event = qnorm(1 - alphah),
    Poisson = qnorm(1 - alphah),
    negbin1 = qnorm(1 - alphah),
    negbin2 = qnorm(1 - alphah))
    
    estVec <- estMat[, 1]
    halfLength <- tailPercentile * estMat[, 2]
    confMat <- matrix(c(estVec -  halfLength, estVec +  halfLength), ncol = 2)    
    
    ## Formatting matrix
    if (formatting)
    {
        colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " "))
        rownames(confMat) <- rownames(estMat)
    }

    return(confMat)    
}
