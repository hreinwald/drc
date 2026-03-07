#' Simulating a dose-response curve
#'
#' Simulation of a dose-response curve with user-specified dose values and error distribution.
#'
#' The distribution for the dose values can either be a fixed set of dose values (a numeric
#' vector) used repeatedly for creating all curves or be a distribution specified as a
#' character string resulting in varying dose values from curve to curve.
#'
#' The error distribution for the response values can be any continuous distribution
#' like \code{\link{rnorm}} or \code{\link{rgamma}}. Alternatively it can be the binomial
#' distribution \code{\link{rbinom}}.
#'
#' @param nosim numeric. The number of simulated curves to be returned.
#' @param fct list. Any built-in function in the package \emph{drc} or a list with similar
#'   components.
#' @param mpar numeric. The model parameters to be supplied to \code{fct}.
#' @param xerror numeric or character. The distribution for the dose values.
#' @param xpar numeric vector supplying the parameter values defining the distribution for the
#'   dose values. If \code{xerror} is a distribution then remember that the number of dose
#'   values also is part of this argument (the first argument).
#' @param yerror numeric or character. The error distribution for the response values.
#' @param ypar numeric vector supplying the parameter values defining the error distribution
#'   for the response values.
#' @param onlyY logical. If TRUE then only the response values are returned (useful in
#'   simulations). Otherwise both dose values and response values (and for binomial data also
#'   the weights) are returned.
#'
#' @return A list with up to 3 components (depending on the value of the \code{onlyY} argument).
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Simulating normally distributed dose-response data
#'
#' ## Model fit to simulate from
#' ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
#'
#' ## 10 random dose-response curves based on the model fit
#' sim10a <- rdrm(10, LL.4(), coef(ryegrass.m1), xerror = ryegrass$conc)
#' sim10a
#'
#' ## Simulating binomial dose-response data
#'
#' ## Model fit to simulate from
#' deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")
#'
#' ## 10 random dose-response curves
#' sim10b <- rdrm(10, LL.2(), coef(deguelin.m1), deguelin$dose, yerror="rbinom", ypar=deguelin$n)
#' sim10b
#'
#' @keywords models nonlinear
"rdrm" <- function(nosim, fct, mpar, xerror, xpar = 1, yerror = "rnorm", ypar = c(0, 1), 
onlyY = FALSE)
{        
    ## Constructing the predictor values
    if (is.numeric(xerror))
    {
        x <- xerror
    } else {
        xFun <- match.fun(xerror)
        x <- do.call(xFun, as.list(xpar))
    }
    lenx <- length(x)
    x <- sort(x)
    x <- rep(x, nosim)
    xMat <- matrix(x, nosim, lenx, byrow = TRUE)
    
    ## Constructing the mean dose-response
    meanVec <- fct$fct(x, matrix(mpar, lenx*nosim, length(mpar), byrow = TRUE))
    
    ## Constructing the simulated response values
    yFun <- match.fun(yerror)
    if (yerror == "rbinom")
    {
        if (length(ypar) == 1)
        {
            ypar <- rep(ypar, lenx*nosim)
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        } else {
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        }
        errorVec <- yFun(lenx*nosim, ypar, meanVec)
        
        yMat <- matrix(errorVec, nosim, lenx, byrow = TRUE)

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {    
            return(list(x = xMat, w = wMat, y = yMat))
        }
    }  else {
        errorVec <- do.call(yFun, c(list(lenx*nosim), as.list(ypar)))
        
        yMat <- matrix(meanVec, nosim, lenx, byrow = TRUE) + errorVec

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {
            return(list(x = xMat, y = yMat))
        }
    }
}
