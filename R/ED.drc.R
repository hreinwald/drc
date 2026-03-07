#' @title Estimating effective doses
#' @description Estimates effective concentration or doses for specified response levels.
#' @param object an object of class 'drc'.
#' @param ... additional arguments passed to methods.
#' @keywords models nonlinear
ED <- function (object, ...) UseMethod("ED", object)

#' @title Estimating effective doses
#'
#' @description
#' \code{ED} estimates effective concentration or doses for one or more specified absolute or
#' relative response levels.
#'
#' @param object an object of class 'drc'.
#' @param respLev a numeric vector containing the response levels.
#' @param interval character string specifying the type of confidence intervals to be supplied.
#'   The default is \code{"none"}. See Details below for more explanation.
#' @param clevel character string specifying the curve id in case estimates for a specific curve
#'   or compound are requested. By default estimates are shown for all curves.
#' @param level numeric. The level for the confidence intervals. The default is 0.95.
#' @param reference character string. Is the upper limit or the control level the reference?
#' @param type character string. Whether the specified response levels are absolute or relative (default).
#' @param lref numeric value specifying the lower limit to serve as reference.
#' @param uref numeric value specifying the upper limit to serve as reference (e.g., 100\%).
#' @param bound logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful
#'   for hormesis models.
#' @param vcov. function providing the variance-covariance matrix or a variance-covariance matrix.
#'   \code{\link{vcov}} is the default, but \code{sandwich} is also an option (for obtaining robust
#'   standard errors).
#' @param display logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).
#' @param logBase numeric. The base of the logarithm in case logarithm transformed dose values are used.
#' @param multcomp logical to switch on output for use with the package \pkg{multcomp} (which needs
#'   to be activated first). Default is FALSE.
#' @param intType string specifying the type of interval to use with the predict method in case the
#'   type of confidence interval chosen is inverse regression.
#' @param ... additional arguments passed to the ED function in the model.
#'
#' @return An invisible matrix containing the estimates and the corresponding estimated standard
#'   errors and possibly lower and upper confidence limits. Or, alternatively, a list with elements
#'   that may be plugged directly into \code{parm} in the package \pkg{multcomp} (when \code{multcomp}
#'   is TRUE).
#'
#' @details
#' There are several options for calculating confidence intervals through the argument \code{interval}.
#' The option \code{"delta"} results in asymptotical Wald-type confidence intervals (using the delta
#' method and the normal or t-distribution depending on the type of response). The option \code{"fls"}
#' produces (possibly skewed) confidence intervals through back-transformation from the logarithm
#' scale (only meaningful in case the parameter in the model is log(ED50) as for the
#' \code{\link{llogistic2}} models). The option \code{"tfls"} is for transforming back and forth from
#' log scale (experimental). The option \code{"inv"} results in confidence intervals obtained through
#' inverse regression.
#'
#' For hormesis models (\code{\link{braincousens}} and \code{\link{cedergreen}}), the additional
#' arguments \code{lower} and \code{upper} may be supplied. These arguments specify the lower and
#' upper limits of the bisection method used to find the ED values.
#'
#' @seealso \code{\link{EDcomp}} for estimating differences and ratios of ED values,
#'   \code{\link{compParm}} for comparing other model parameters, and \code{\link{backfit}}.
#'
#' @examples
#' ## Fitting 4-parameter log-logistic model
#' ryegrass.m1 <- drm(ryegrass, fct = LL.4())
#'
#' ## Calculating EC/ED values
#' ED(ryegrass.m1, c(10, 50, 90))
#'
#' ## Also displaying 95% confidence intervals
#' ED(ryegrass.m1, c(10, 50, 90), interval = "delta")
#'
#' @author Christian Ritz
#' @keywords models nonlinear
"ED.drc" <- function(object, 
                     respLev, 
                     interval = c("none", "delta", "fls", "tfls", "inv"), 
                     clevel = NULL,
                     level = ifelse(!(interval == "none"), 0.95, NULL), 
                     reference = c("control", "upper"), 
                     type = c("relative", "absolute"), 
                     lref, uref, bound = TRUE, 
                     vcov. = vcov,
                     display = TRUE, 
                     logBase = NULL, 
                     multcomp = FALSE, 
                     intType = "confidence", ...)
{
    interval <- match.arg(interval)
    reference <- match.arg(reference)
    type <- match.arg(type)
    
    ## Checking 'respLev' vector: it should be numbers between 0 and 100
    if ( (type == "relative") && (bound) ) 
    {
        if (any(respLev <= 0 | respLev >= 100)) 
        {
            stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
        }
    }

    ## Retrieving relevant quantities
    EDlist <- object$fct[["edfct"]]  
    if (is.null(EDlist)) {stop("ED values cannot be calculated")}         
    indexMat <- object[["indexMat"]]
    parmMat <- object[["parmMat"]]

    curveNames <- colnames(parmMat)  # colnames(object$"parmMat")
    options(warn = -1)  # switching off warnings caused by coercion in the if statement
    if (any(is.na(as.numeric(curveNames))))
    {
        curveOrder <- order(curveNames)
    } else { # if names are numbers then skip re-ordering
        curveOrder <- 1:length(curveNames)
    }
    options(warn = 0)  # normalizing behaviour of warnings
    
    strParm0 <- curveNames[curveOrder]
    indexMat <- indexMat[, curveOrder, drop = FALSE]
    parmMat <- parmMat[, curveOrder, drop = FALSE]
    
    strParm <- strParm0
    #vcMat <- vcov.(object)
    if (is.function(vcov.))  # following a suggestion by Andrea Onofri
    {
      vcMat <- vcov.(object)
    } else {
      vcMat <- vcov.
    }
    
    ## Defining vectors and matrices needed
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM    
#    lenEB <- ncolIM    
    lenPV <- length(respLev)  # used twice below
    noRows <- ncolIM * lenPV
    dimNames <- rep("", noRows)  # lenEB*lenPV, 2)
    EDmat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)
    oriMat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)

    ## Skipping curve id if only one curve is present
    if (identical(length(unique(strParm)), 1)) 
    {
        strParm[indexVec] <- rep("", ncolIM)
    } else {
        strParm <- paste(strParm, ":", sep = "")
    }

    ## Calculating estimates and estimated standard errors
    rowIndex <- 1
    lenIV <- length(indexVec)
    dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
    intMat <- NULL
    for (i in indexVec)
    {
        parmChosen <- parmMat[, i]
        parmInd <- indexMat[, i]
        varCov <- vcMat[parmInd, parmInd]

        if ((is.null(clevel)) || (strParm0[i] %in% clevel))
        {
        for (j in 1:lenPV)
        {
            EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)            
            EDval <- EDeval[[1]]
            dEDval <- EDeval[[2]]
            dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval 
            
            oriMat[rowIndex, 1] <- EDval
            oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
                   
            if (!is.null(logBase))
            {
                EDval <- logBase^(EDval)                
                dEDval <- EDval * log(logBase) * dEDval
            }
            EDmat[rowIndex, 1] <- EDval
            EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)

            dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
            rowIndex <- rowIndex + 1
        }
            if (interval == "inv")
            {
                intMat <- rbind(intMat, t(EDinvreg1(object, respLev, strParm0[i], 
                                                    intType = intType, level = level, type = type)))
            }
          
        } else {
            rowsToRemove <- rowIndex:(rowIndex + lenPV - 1)
            EDmat <- EDmat[-rowsToRemove, , drop = FALSE]
            dimNames <- dimNames[-rowsToRemove]
        }
        
    }
    
    ## Defining column names
    colNames <- c("Estimate", "Std. Error")
    
    ## Calculating the confidence intervals
    if (interval == "delta")
    {
        intMat <- confint.basic(EDmat, level, object$"type", df.residual(object), FALSE)
        intLabel <- "Delta method"
    }
    
    if (interval == "tfls")
    {
        intMat <- exp(confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2] / oriMat[, 1]), ncol = 2), 
        level, object$"type", df.residual(object), FALSE))
        intLabel <- "To and from log scale"       
    }

    if (interval == "fls")
    {        
        if (is.null(logBase)) 
        {
            logBase <- exp(1)
            EDmat[, 1] <- exp(EDmat[, 1])  # back-transforming log ED values
        }

        intMat <- logBase^(confint.basic(oriMat, level, object$"type", df.residual(object), FALSE))
        intLabel <- "Back-transformed from log scale"  

        ## Dropping estimated standard errors (not relevant after back transformation)        
        EDmat <- EDmat[, -2, drop = FALSE]  
        colNames <- colNames[-2]
#        colNames <- c(colNames[-2], "Lower", "Upper")  # standard errors not relevant        
    }
    
    if (interval == "inv")
    {
        EDmat <- EDmat[, -2, drop = FALSE]  
        colNames <- colNames[-2]
        intLabel <- "Inverse regression"  
    }
    
    if (identical(interval, "none"))
    {
        intLabel <- NULL
    } else {
        EDmat <- as.matrix(cbind(EDmat, intMat))
        colNames <- c(colNames, "Lower", "Upper")         
    } 
    dimnames(EDmat) <- list(dimNames, colNames)
    rownames(EDmat) <- paste("e", rownames(EDmat), sep = ":")
    resPrint(EDmat, "Estimated effective doses", interval, intLabel, display = display)
    
    if(multcomp)
    {  
        EDmat1 <- EDmat[, 1]
        namesVec <- names(EDmat1)  # paste("e", names(EDmat1), sep = ":")
#        names(EDmat1) <- namesVec
        
        EDmat1VC <- (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat), drop = FALSE]
        colnames(EDmat1VC) <- namesVec
        rownames(EDmat1VC) <- namesVec
        
        invisible(list(#EDdisplay = EDmat, 
#                       EDmultcomp = parm(EDmat[, 1], (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat), drop = FALSE])))
                        EDmultcomp = parm(EDmat1, EDmat1VC)))
    } else {
        invisible(EDmat)     
    }   
}

