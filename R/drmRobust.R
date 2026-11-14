#' @title Robust estimation functions for drm
#' @param robust character string, the rho function (see \code{\link{drm}}).
#' @param fctCall the (matched) call of \code{drm()}.
#' @param lenData number of observations.
#' @param lenPar number of estimated parameters.
#' @param callEnv environment in which \code{fctCall} is evaluated for the
#'   auxiliary LTS fit (scale estimate of \code{"trimmed"}, \code{"tukey"} and
#'   \code{"winsor"}); \code{drm()} passes the environment it was called from.
#' @keywords internal
"drmRobust" <- function(robust, fctCall, lenData, lenPar, callEnv = parent.frame())
{ 
    validRobust <- c("mean", "median", "trimmed", "tukey", "winsor", "lms", "lts")
    if (!is.character(robust) || length(robust) != 1L || !(robust %in% validRobust))
    {
        stop("Invalid 'robust' argument: must be one of ", paste0("'", validRobust, "'", collapse = ", "), call. = FALSE)
    }

    ## Finding robust scale estimate for trimmed and winsorised means and tukey
    if (robust%in%c("trimmed", "tukey", "winsor"))
    {
        call1 <- fctCall
        call1$"robust" <- "lts"

        scaleEst <- mad(residuals(eval(call1, callEnv)), 0)
    }
    

    ## Defining distance functions
    quadratic <- function(x) {x*x}
    
    lms <- function(x) {median(x*x)}

    noRes <- floor((lenData+lenPar+1)/2)
    lts <- function(x) {sum(((x[order(x)])[1:noRes])^2)}
 
    metricTrim <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        cVal <- 1.345
        retVal <- x*x
        retVal[abs(x) > cVal] <- cVal^2        
        retVal
    }

    metricWinsor <- function(x) 
    {       
        if (all(is.na(x))) {return(x)}
    
        cVal <- 1.345

        x <- x/scaleEst

        retVal <- x*x

        indexVec <- !is.na(x) & abs(x) > cVal  # NA-safe subsetting (after upstream, thanks to Gary Carnahan)
        retVal[indexVec] <- (cVal * (2 * abs(x) - cVal))[indexVec]
        retVal
    }

    tukeyBiweight <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        Rval <- 4.685
        retVal <- (x^6)/(Rval^4) - 3*(x^4)/(Rval^2) + 3*x*x
        retVal[abs(x) > Rval] <- Rval * Rval
        retVal
    }


    ## Assigning objective function
    robustFct <- switch(robust, 
                        mean = quadratic, 
                        median = abs, 
                        trimmed = metricTrim, 
                        tukey = tukeyBiweight, 
                        winsor = metricWinsor, 
                        lms = lms, 
                        lts = lts)
    
    return(robustFct)    
}
