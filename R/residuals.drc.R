#' @title Extracting residuals from the fitted dose-response model
#'
#' @description
#' \code{residuals} extracts different types of residuals from an object of
#' class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param typeRes character string specifying the type of residual to be
#'   returned: raw/working residuals, residuals standardised using the
#'   estimated residual standard error, or studentised residuals based on the
#'   H matrix of partial derivatives of the model function.
#' @param trScale logical value indicating whether or not to return residuals
#'   on the transformed scale (in case a Box-Cox transformation was applied).
#' @param ... additional arguments.
#'
#' @return The raw (also called working) residuals or some kind of scaled
#'   residuals extracted from \code{object}.
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' ## Displaying the residual plot (raw residuals)
#' plot(fitted(ryegrass.m1), residuals(ryegrass.m1))
#'
#' ## Using the standardised residuals
#' plot(fitted(ryegrass.m1), residuals(ryegrass.m1, typeRes = "standard"))
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"residuals.drc" <-
function(object, typeRes = c("working", "standardised", "studentised"), trScale = TRUE, ...)
{
    typeRes <- match.arg(typeRes)
    
#    rawResiduals <- object$"dataList"$"resp" - fitted(object)  # changed 29/12 2012

    if (trScale && (!is.null(object$"boxcox")))
    {
        ## Defining Box-Cox transformation function
        bcfct <- function(x, lambda, bctol = 0.02, add)
        # same function as in drm()
        {
            if (abs(lambda) > bctol)
            {
                return(((x + add)^lambda - 1)/lambda)
            } else {
                return(log(x + add))    
            }
        }
        objBC <- object$"boxcox"
        rawResiduals <- object$"dataList"$"resp" - bcfct(fitted(object), lambda = objBC$"lambda", add = objBC$"bcAdd")
    } else {
        rawResiduals <- object$"dataList"$"origResp" - fitted(object)
    }
    
    if (identical(typeRes, "standardised"))
    {
        rstan <- object$"estMethod"$"rstanfct"
        if (is.null(rstan)) 
        {
            cat("Scale parameter fixed at 1. So working residuals are returned\n\n")
#            return(object$"predres"[, 2])
            return(rawResiduals)
            
#            stop("No standardisation available")
        } else {
#            return(object$"predres"[, 2] / rstan(object))  
            return(rawResiduals / rstan(object))  
            
#            return( object$"predres"[, 2] / sqrt(summary(object)$"resVar") )
        }
    } 
    
    if (identical(typeRes, "studentised"))
    {
        Xmat <- object$"deriv1"
        if (is.null(Xmat)) {stop("Studentised residuals not available")}
        
        Xprod <- solve(t(Xmat) %*% Xmat)
        nrowX <- nrow(Xmat)        
        Hdiag <- rep(NA, nrowX)
        for (i in 1:nrowX)
        {
            Hdiag[i] <- Xmat[i, ] %*% Xprod %*% t(Xmat[i, , drop = FALSE])
        }
#        print(length(Hdiag))
#        print(dim(object$"predres"))
        scaleEst0 <- summary(object)$"resVar"
        scaleEst <- ifelse(is.na(scaleEst0), 1, scaleEst0)  
        # to handle response types that are not continuous/quantitative
        
#        return(object$"predres"[, 2] / sqrt(scaleEst * (1 - Hdiag)))
        return(rawResiduals / sqrt(scaleEst * (1 - Hdiag)))
    }
    
    if (identical(typeRes, "working"))
    {
#        return(object$"predres"[, 2])
        return(rawResiduals)
    }
}
