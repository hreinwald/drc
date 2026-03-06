#' Five-Parameter Log-Logistic Model with log(ED50) as Parameter
#'
#' A five-parameter log-logistic model where the ED50 is parameterised on the
#' log scale. The mean function is:
#' \deqn{f(x) = c + \frac{d - c}{(1 + \exp(b(\log(x) - e)))^f}}{f(x) = c + (d-c)/(1+exp(b(log(x)-e)))^f}
#' where \code{e} is the logarithm of the ED50 (not exponentiated).
#'
#' @param fixed numeric vector of length 5. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 5 giving the names of the
#'   parameters \code{b}, \code{c}, \code{d}, \code{e}, and \code{f}.
#' @param ss character string indicating the self-starter version to use.
#'   One of \code{"1"} (default), \code{"2"}, or \code{"3"}.
#' @param ssfct optional self-starter function. If provided, overrides the
#'   built-in self-starter selected by \code{ss}.
#' @param fctName optional character string specifying the name of the function.
#' @param fctText optional character string providing a short description of
#'   the function.
#'
#' @return A list of class \code{"llogistic"} containing the nonlinear function,
#'   self-starter function, parameter names, and related helper functions.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{llogistic}}, \code{\link{LL2.2}}, \code{\link{LL2.3}},
#'   \code{\link{LL2.4}}, \code{\link{LL2.5}}
#'
#' @keywords models nonlinear
"llogistic2" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ss = c("1", "2", "3"), ssfct = NULL,
fctName, fctText)
{
    ## Matching 'adjust' argument
    ss <- match.arg(ss)
    
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the basic non-linear function
    bfct <- function(x, parm)
    {
        parm[2] + (parm[3]-parm[2])/((1+(x/exp(parm[4]))^parm[1]))^parm[5]
    }

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose)-parmMat[, 4])))^parmMat[, 5])
    }


    ## Defining the self starter function
    
    ## Defining self starter function based on supplied function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        
    ## Version 1 (default)    
    if (ss == "1")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

            startVal <- rep(0, numParm)

            startVal[3] <- max(y) + 0.001  # the d parameter
            startVal[2] <- min(y) - 0.001  # the c parameter
            startVal[5] <- 1  # better choice may be possible!        
        
            if (length(unique(x))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 

            indexT2 <- (x > 0 & is.finite(x))
#            if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
            x2 <- x[indexT2]
            y2 <- y[indexT2]

            startVal[c(1,4)] <- find.be2(x2, y2, startVal[2] - 0.001, startVal[3])
            startVal[4] <- log(startVal[4])
            # 0.001 to avoid 0 in the denominator

            return(startVal[notFixed])
        }
    }

    ## Version 2
    if (ss == "2")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

            dVal <- ifelse(notFixed[3], 1.01*max(y), fixed[3])        
            cVal <- ifelse(notFixed[2], 0.99*min(y), fixed[2])

            fVal <- 1  # need not be updated with value in 'fixed[5]'
            # better choice than 1 may be possible!        
        
            if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 

            # Cutting away response values close to d
            indexT1a <- x > 0
            x2 <- x[indexT1a]
            y2 <- y[indexT1a]
            
            beVec <- find.be2(x2, y2, cVal, dVal)
            bVal <- beVec[1]
            eVal <- log(beVec[2])
            
            return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        }
    }

    ## Version 3
    if (ss == "3")
    {
        ssfct <- function(dframe)
        {
            x <- dframe[, 1]
            y <- dframe[, 2]

            cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2])
            dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3])
            fVal <- 1  # need not be updated with value in 'fixed[5]'
        
            if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
            # only estimate of upper limit if a single unique dose value 
           
            beVec <- find.be1(x, y, cVal, dVal)
            bVal <- beVec[1]
            eVal <- log(beVec[2])
            
            return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        }
    }
    }

    ## Finding b and e based on stepwise increments
    find.be1 <- function(x, y, c, d)
    {
        unix <- unique(x)
        uniy <- tapply(y, x, mean)
        lenx <- length(unix)
        
        j <- 2
        for (i in 2:lenx)
        {
            crit1 <- (uniy[i] > (d + c)/2) && (uniy[i-1] < (d + c)/2)
            crit2 <- (uniy[i] < (d + c)/2) && (uniy[i-1] > (d + c)/2)
            if (crit1 || crit2) break
            j <- j + 1
        }
        eVal <- (unix[j] + unix[j-1])/2
        bVal <- -sign(uniy[j] - uniy[j-1])  # -(uniy[j] - uniy[j-1]) / (unix[j] - unix[j-1])
        return(as.vector(c(bVal, eVal)))  
    }
    
    ## Finding b and e based on linear regression
    find.be2 <- function(x, y, c, d)
    {
        logitTrans <- log((d - y)/(y - c))  

        lmFit <- lm(logitTrans ~ log(x))

        coefVec <- coef(lmFit)
        bVal <- coefVec[2]        
        eVal <- exp(-coefVec[1]/bVal)    

        return(as.vector(c(bVal, eVal)))
    }

   
    ## Defining names
    names <- names[notFixed]
    
    
    ##Defining the first derivatives (in the parameters)        
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(log(dose) - parmMat[, 4]))
        t5 <- (1 + t2)^parmMat[, 5]                  

        cbind( -t1 * xlogx(dose/exp(parmMat[, 4]), parmMat[, 1], parmMat[, 5] + 1) * parmMat[, 5], 
               1 - 1/t5, 
               1/t5, 
               t1 * parmMat[, 5] * divAtInf(t2, (1 + t2)^(parmMat[, 5] + 1)) * parmMat[, 1], 
               -t1 * divAtInf(log(1+t2), t5) )[, notFixed]
    }
    deriv2 <- NULL


    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        temp1 <- exp(parmMat[, 1]*(log(x) - parmMat[, 4]))  # x/parmMat[, 4]          
        temp2 <- 1 + temp1
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*temp1*parmMat[, 1]/x
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4 
    }


#    ## Setting the limits
#    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    ## (returning ED values and corresponding standard errors on log scale)
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
#        if (type == "absolute") 
#        {
#            p <- 100*((parmVec[3] - respl)/(parmVec[3] - parmVec[2]))
#        } else {  
#            p <- respl
#        }
#        if ( (parmVec[1] < 0) && (reference == "control") )
#        {
#            p <- 100 - p
#        }
        p <- EDhelper(parmVec, respl, reference, type)        
    
        tempVal1 <- 100/(100-p)
        tempVal2 <- log(tempVal1^(1/parmVec[5]) - 1)
        lEDp <- parmVec[4] + tempVal2 / parmVec[1]

        lEDder <- 
        c(-tempVal2/(parmVec[1]^2), 
        0, 0, 1, 
        tempVal1^(1/parmVec[5]-1)/(parmVec[1]*parmVec[5]*(tempVal1^(1/parmVec[5]-1))))

        return(list(lEDp, lEDder[notFixed]))
    }
 
    
    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + parmVec[4])
    } 
    
    ## Defining functions returning lower and upper limit and monotonicity
    lowerAs <- pickParm(parmVec, notFixed, 2)
    upperAs <- pickParm(parmVec, notFixed, 3)
    monoton <- monoParm(parmVec, notFixed, 1, -1)    
    
    ## Setting function details
    if (missing(fctName))
    {
        fctName <- as.character(match.call()[[1]])
    }  
    if (missing(fctText))
    {
        fctText <- "Log-logistic (log(ED50) as parameter)"
    }         
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, bfct = bfct, inversion = invfct,
#    lowerc=lowerLimits, upperc=upperLimits,
    name = fctName, 
    text = fctText, 
    noParm = sum(is.na(fixed)),
    lowerAs = lowerAs, 
    upperAs = upperAs, 
    monoton = monoton)
    
    class(returnList) <- "llogistic"
    invisible(returnList)
}

#' Construct Text for Model with Fixed Lower and Upper Limits
#'
#' Helper function that appends lower and upper limit information to a model
#' description string.
#'
#' @param modelStr character string with the base model description.
#' @param upper numeric value for the fixed upper limit.
#'
#' @return A character string describing the model with its fixed limits.
#'
#' @keywords internal
lowupFixed <- function(modelStr, upper)
{
    paste(modelStr, "with lower limit at 0 and upper limit at", upper)
}

#' Construct Text for Model with Fixed Lower Limit
#'
#' Helper function that appends lower limit information to a model description
#' string.
#'
#' @param modelStr character string with the base model description.
#'
#' @return A character string describing the model with its fixed lower limit.
#'
#' @keywords internal
lowFixed <- function(modelStr)
{
    paste(modelStr, "with lower limit at 0")
}

#' Construct Text for Model with Fixed Upper Limit
#'
#' Helper function that appends upper limit information to a model description
#' string.
#'
#' @param modelStr character string with the base model description.
#' @param upper numeric value for the fixed upper limit.
#'
#' @return A character string describing the model with its fixed upper limit.
#'
#' @keywords internal
upFixed <- function(modelStr, upper)
{
    paste(modelStr, "with upper limit at", upper)
}

#' Two-Parameter Log-Logistic Model with log(ED50) as Parameter
#'
#' A two-parameter log-logistic model where the lower limit is fixed at 0 and
#' the upper limit is fixed at a specified value (default 1). The estimated
#' parameters are the slope \code{b} and the log(ED50) \code{e}.
#'
#' @param upper numeric value giving the fixed upper limit. Defaults to 1.
#' @param fixed numeric vector of length 2. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 2 giving the names of the
#'   parameters \code{b} and \code{e}.
#' @param \dots additional arguments passed to \code{\link{llogistic2}}.
#'
#' @return A list of class \code{"llogistic"} with the nonlinear function,
#'   self-starter, and related components.
#'
#' @seealso \code{\link{llogistic2}}, \code{\link{LL2.3}}, \code{\link{LL2.4}},
#'   \code{\link{LL2.5}}
#'
#' @examples
#' earthworms.m1 <- drm(number/total ~ dose, weights = total,
#'   data = earthworms, fct = LL2.2(), type = "binomial")
#'
#' @keywords models nonlinear
"LL2.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1], 0, upper, fixed[2], 1), 
    names = c(names[1], "c", "d", names[2], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Log-logistic (log(ED50) as parameter)", upper),   
    ...) )
}

#' Three-Parameter Log-Logistic Model with log(ED50) and Lower Limit at 0
#'
#' A three-parameter log-logistic model where the lower limit is fixed at 0.
#' The estimated parameters are the slope \code{b}, the upper limit \code{d},
#' and the log(ED50) \code{e}.
#'
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 3 giving the names of the
#'   parameters \code{b}, \code{d}, and \code{e}.
#' @param \dots additional arguments passed to \code{\link{llogistic2}}.
#'
#' @return A list of class \code{"llogistic"} with the nonlinear function,
#'   self-starter, and related components.
#'
#' @seealso \code{\link{llogistic2}}, \code{\link{LL2.2}}, \code{\link{LL2.4}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.3())
#'
#' @keywords models nonlinear
"LL2.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1], 0, fixed[2:3], 1), 
    names = c(names[1], "c", names[2:3], "f"),
    fctName = as.character(match.call()[[1]]),
    fctText = lowFixed("Log-logistic (log(ED50) as parameter)"), 
    ...) )
}

#' Three-Parameter Log-Logistic Model with log(ED50) and Fixed Upper Limit
#'
#' A three-parameter log-logistic model where the upper limit is fixed at a
#' specified value (default 1). The estimated parameters are the slope \code{b},
#' the lower limit \code{c}, and the log(ED50) \code{e}.
#'
#' @param upper numeric value giving the fixed upper limit. Defaults to 1.
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 3 giving the names of the
#'   parameters \code{b}, \code{c}, and \code{e}.
#' @param \dots additional arguments passed to \code{\link{llogistic2}}.
#'
#' @return A list of class \code{"llogistic"} with the nonlinear function,
#'   self-starter, and related components.
#'
#' @seealso \code{\link{llogistic2}}, \code{\link{LL2.2}}, \code{\link{LL2.3}}
#'
#' @keywords models nonlinear
"LL2.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed[1:2], upper, fixed[3], 1), 
    names = c(names[1:2], "d", names[3], "f"), 
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Log-logistic (log(ED50) as parameter)", upper), 
    ...) )
}

#' Four-Parameter Log-Logistic Model with log(ED50) as Parameter
#'
#' A four-parameter log-logistic model where the ED50 is parameterised on the
#' log scale. The asymmetry parameter \code{f} is fixed at 1. The estimated
#' parameters are the slope \code{b}, the lower limit \code{c}, the upper
#' limit \code{d}, and the log(ED50) \code{e}.
#'
#' @param fixed numeric vector of length 4. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 4 giving the names of the
#'   parameters \code{b}, \code{c}, \code{d}, and \code{e}.
#' @param \dots additional arguments passed to \code{\link{llogistic2}}.
#'
#' @return A list of class \code{"llogistic"} with the nonlinear function,
#'   self-starter, and related components.
#'
#' @seealso \code{\link{llogistic2}}, \code{\link{LL2.3}}, \code{\link{LL2.5}}
#'
#' @examples
#' spinach.m1 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL2.4())
#'
#' @keywords models nonlinear
"LL2.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic2(fixed = c(fixed, 1), names = c(names, "f"),
    fctName = as.character(match.call()[[1]]), ...) )
}

#' Five-Parameter Generalised Log-Logistic Model with log(ED50) as Parameter
#'
#' A five-parameter generalised log-logistic model where the ED50 is
#' parameterised on the log scale. All five parameters (\code{b}, \code{c},
#' \code{d}, \code{e}, \code{f}) are estimated.
#'
#' @param fixed numeric vector of length 5. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated.
#' @param names character vector of length 5 giving the names of the
#'   parameters \code{b}, \code{c}, \code{d}, \code{e}, and \code{f}.
#' @param \dots additional arguments passed to \code{\link{llogistic2}}.
#'
#' @return A list of class \code{"llogistic"} with the nonlinear function,
#'   self-starter, and related components.
#'
#' @seealso \code{\link{llogistic2}}, \code{\link{LL2.3}}, \code{\link{LL2.4}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.5())
#'
#' @keywords models nonlinear
"LL2.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return( llogistic2(fixed = fixed, names = names, 
    fctName = as.character(match.call()[[1]]),
    fctText = "Generalised log-logistic (log(ED50) as parameter)", ...) )
}

