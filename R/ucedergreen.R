#' @title U-shaped Cedergreen-Ritz-Streibig model
#'
#' @description
#' \code{ucedergreen} provides a very general way of specifying the Cedergreen-Ritz-Streibig
#' modified log-logistic model for describing u-shaped hormesis, under various constraints on the parameters.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed.
#'   NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters (should not contain ":").
#'   The order of the parameters is: b, c, d, e, f.
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used.
#' @param alpha numeric value between 0 and 1, reflecting the steepness of the hormesis peak.
#'   This argument must be specified.
#'
#' @details
#' The u-shaped model is given by the expression
#' \deqn{f(x) = c + d - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' @return A list containing the non-linear function, the self starter function
#'   and the parameter names.
#'
#' @references
#' Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005)
#' Improved empirical models describing hormesis,
#' \emph{Environmental Toxicology and Chemistry} \bold{24}, 3166--3172.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{UCRS.4a}}, \code{\link{UCRS.5a}}, \code{\link{drm}}
#'
#' @keywords models nonlinear
"ucedergreen" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
alpha)
{
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}
    
    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec

    ## Defining the function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        numTerm <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*exp(-1/dose^alpha)
        denTerm <- 1 + exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        parmMat[, 3] - numTerm/denTerm
    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)   
            initval[1] <- -initval[1]
            initval[5] <- 0  # better solution?
    
            return(initval[notFixed])
        }
    }   
    
#    ## Setting the names of the parameters
#    names <- names[notFixed]


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }


    ## Defining derivatives
    
#    ## Constructing a helper function
#    xlogx <- function(x, p)
#    {
#        lv <- (x < 1e-12)
#        nlv <- !lv
#        
#        rv <- rep(0, length(x))
#        
#        xlv <- x[lv] 
#        rv[lv] <- log(xlv^(xlv^p[lv]))
#        
#        xnlv <- x[nlv]
#        rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
#    
#        rv
#    }
    
    ## Specifying the derivatives    
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t0 <- exp(-1/(dose^alpha))
        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*t0
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1/t3, 
               1 - 1/t3, 
               -t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               -t0/t3 )[, notFixed]
    }
        
    deriv2 <- NULL

    ## Setting limits
#    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

    ## Defining the ED function
    edfct <- function(parm, p, lower = 1e-4, upper = 10000, ...)
    {    
        cedergreen(fixed =  fixed, names = names, alpha = alpha)$edfct(parm, 100 - p, lower, upper, ...) 
    }

#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair, upper = 10000, interval = c(1e-4, 10000))
#    {
#        cedergreen(alpha = alpha)$sifct(parm1, parm2, 100-pair, upper, interval)
#    }    

    ## Finding the maximal hormesis
    maxfct <- function(parm, upper, interval)
    {
       retVal <- cedergreen(fixed =  fixed, names = names, alpha = alpha)$maxfct(parm, upper, interval)
#       retVal[2] <- (parm[2] + parm[3]) - (retVal[2] - parm[2])
       retVal[2] <- (parm[2] + parm[3]) - retVal[2]
              
       return(retVal)
    }

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], edfct = edfct, maxfct = maxfct,
    name = "ucedergreen",
    text = "U-shaped Cedergreen-Ritz-Streibig", 
    noParm = sum(is.na(fixed)))
    
    class(returnList) <- "UCRS"
    invisible(returnList)
}


#' @title U-shaped CRS model with lower limit 0 (alpha=1)
#'
#' @description
#' Four-parameter u-shaped CRS hormesis model with lower limit fixed at 0 and alpha=1.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.5a}}, \code{\link{CRS.4a}}
#'
#' @keywords models nonlinear
"UCRS.4a" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 1, ...))
}

#' @title Alias for UCRS.4a
#' @description \code{uml3a} is an alias for \code{\link{UCRS.4a}}.
#' @seealso \code{\link{UCRS.4a}}
#' @keywords models nonlinear
uml3a <- UCRS.4a

#' @title U-shaped CRS model with lower limit 0 (alpha=0.5)
#'
#' @description
#' Four-parameter u-shaped CRS hormesis model with lower limit fixed at 0 and alpha=0.5.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.4a}}, \code{\link{CRS.4b}}
#'
#' @keywords models nonlinear
"UCRS.4b" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.5, ...))
}

#' @title Alias for UCRS.4b
#' @description \code{uml3b} is an alias for \code{\link{UCRS.4b}}.
#' @seealso \code{\link{UCRS.4b}}
#' @keywords models nonlinear
uml3b <- UCRS.4b

#' @title U-shaped CRS model with lower limit 0 (alpha=0.25)
#'
#' @description
#' Four-parameter u-shaped CRS hormesis model with lower limit fixed at 0 and alpha=0.25.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.4a}}, \code{\link{CRS.4c}}
#'
#' @keywords models nonlinear
"UCRS.4c" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.25, ...))
}

#' @title Alias for UCRS.4c
#' @description \code{uml3c} is an alias for \code{\link{UCRS.4c}}.
#' @seealso \code{\link{UCRS.4c}}
#' @keywords models nonlinear
uml3c <- UCRS.4c

#' @title U-shaped CRS five-parameter model (alpha=1)
#'
#' @description
#' Five-parameter u-shaped CRS hormesis model with alpha=1.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.4a}}, \code{\link{CRS.5a}}
#'
#' @keywords models nonlinear
"UCRS.5a" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 1, ...))
}

#' @title Alias for UCRS.5a
#' @description \code{uml4a} is an alias for \code{\link{UCRS.5a}}.
#' @seealso \code{\link{UCRS.5a}}
#' @keywords models nonlinear
uml4a <- UCRS.5a

#' @title U-shaped CRS five-parameter model (alpha=0.5)
#'
#' @description
#' Five-parameter u-shaped CRS hormesis model with alpha=0.5.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.5a}}, \code{\link{CRS.5b}}
#'
#' @keywords models nonlinear
"UCRS.5b" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.5, ...))
}

#' @title Alias for UCRS.5b
#' @description \code{uml4b} is an alias for \code{\link{UCRS.5b}}.
#' @seealso \code{\link{UCRS.5b}}
#' @keywords models nonlinear
uml4b <- UCRS.5b

#' @title U-shaped CRS five-parameter model (alpha=0.25)
#'
#' @description
#' Five-parameter u-shaped CRS hormesis model with alpha=0.25.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{ucedergreen}}.
#'
#' @return A list (see \code{\link{ucedergreen}}).
#'
#' @seealso \code{\link{ucedergreen}}, \code{\link{UCRS.5a}}, \code{\link{CRS.5c}}
#'
#' @keywords models nonlinear
"UCRS.5c" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.25, ...))
}

#' @title Alias for UCRS.5c
#' @description \code{uml4c} is an alias for \code{\link{UCRS.5c}}.
#' @seealso \code{\link{UCRS.5c}}
#' @keywords models nonlinear
uml4c <- UCRS.5c
