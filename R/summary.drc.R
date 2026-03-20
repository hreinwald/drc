#' @title Summarising non-linear model fits
#'
#' @description
#' \code{summary} compiles a comprehensive summary for objects of class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param od logical. If TRUE adjustment for over-dispersion is used.
#' @param pool logical. If TRUE curves are pooled. Otherwise they are not. This
#'   argument only works for models with independently fitted curves as
#'   specified in \code{\link{drm}}.
#' @param ... additional arguments.
#'
#' @return A list of summary statistics that includes parameter estimates and
#'   estimated standard errors.
#'
#' @seealso \code{\link{drm}}, \code{\link{coef.drc}}, \code{\link{confint.drc}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' summary(ryegrass.m1)
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"summary.drc" <-
function(object, od = FALSE, pool = TRUE, ...)
{
    ## Calculating variance-covariance matrix from Hessian
    parVec <- as.vector(coef(object))
    varMat <- vcov(object, od = od, pool = pool)
        
    ## Calculating estimated residual variance 
    ## and unscaled variance-covariance matrix
    resVar <- rse(object, TRUE)
    if (!is.na(resVar))
    {
        varMat.us <- varMat / (2*resVar)    
    } else {
        varMat.us <- NULL
    }
   
    ## Calculating the residual standard error(s)
    if ((!is.null(object$"objList")) && (!pool))
    {
        objList <- object$"objList"
        lenol <- length(objList)
        
        rseMat <- matrix(NA, lenol, 2)
        rownames(rseMat) <- names(objList) 
        resVarVec <- as.vector(unlist(lapply(objList, rse, resvar = TRUE)))
        rseMat[, 1] <- sqrt(resVarVec)
        rseMat[, 2] <- as.vector(unlist(lapply(objList, df.residual)))
    } else {
        resVar <- rse(object, TRUE)
        
        rseMat <- matrix(NA, 1, 2)
        rownames(rseMat) <- ""
        rseMat[1, 1] <- sqrt(resVar)
        rseMat[1, 2] <- df.residual(object)
    }
    colnames(rseMat) <- c("rse", "df")
    
    diagVar <- diag(varMat)
    estSE <- numeric(length(diagVar))
    validVar <- diagVar >= 0
    estSE[validVar] <- sqrt(diagVar[validVar])
    estSE[!validVar] <- NaN

    ## Calculating estimated standard errors for robust methods
    
    ## M-estimators
    if (!is.null(object$robust) && object$robust%in%c("metric trimming", "metric Winsorizing", "Tukey's biweight"))
    {
         # Observed "information"-type of variance-covariance matrix
         estSE <- sqrt(resVar * diag(solve(object[["fit"]][["hessian"]])))
    }


    ## Forming a matrix of results        
    parNames <- object$"parNames"[[1]]    
    resultMat <- matrix(NA, length(parVec), 4, 
    dimnames = list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))    

    resultMat[, 1] <- parVec
    resultMat[, 2] <- estSE
    tempStat <- resultMat[, 1] / resultMat[, 2]
    resultMat[, 3] <- tempStat
    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (object$"type" == "continuous")
    {
        pFct <- function(x) {pt(x, df.residual(object))}
    } else {
        pFct <- pnorm
    }    
    resultMat[, 4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))

    fctName <- deparse(object$call$fct)    

    sumObj <- list(resVar, varMat, resultMat, object$"boxcox", fctName, object$"robust", NULL, object$"type", 
    df.residual(object), varMat.us, object$"fct"$"text", object$"fct"$"noParm", rseMat)
    names(sumObj) <- c("resVar", "varMat", "coefficients", "boxcox", "fctName", "robust", "varParm", "type", 
    "df.residual", "cov.unscaled", "text", "noParm", "rseMat")
    
    class(sumObj) <- c("summary.drc")
    return(sumObj)
}
