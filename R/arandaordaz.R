#' Asymptotic Regression Model
#'
#' The base function for the asymptotic regression model, providing the mean
#' function and self starter for a three-parameter model.
#'
#' The asymptotic regression model is a three-parameter model with mean function:
#'
#' \deqn{f(x) = c + (d-c)(1-\exp(-x/e))}
#'
#' The parameter \eqn{c} is the lower limit (at \eqn{x=0}), \eqn{d} is the upper limit,
#' and \eqn{e>0} determines the steepness of the increase.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names vector of character strings giving the names of the parameters
#'   (should not contain ":").
#' @param fctName optional character string used internally by convenience functions.
#' @param fctText optional character string used internally by convenience functions.
#'
#' @return A list of class \code{drcMean}, containing the mean function, the self starter
#'   function, the parameter names, and other components such as derivatives and a function
#'   for calculating ED values.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{AR.2}}, \code{\link{AR.3}}, \code{\link{EXD.2}}, \code{\link{EXD.3}}
#'
#' @keywords models nonlinear
"arandaordaz" <- function(
fixed = c(NA, NA, NA), names = c("a", "b", "c"), fctName, fctText)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[, 1] + (parmMat[, 2] - parmMat[, 1]) * (1 - exp( -parmMat[, 3] * dose))
    }

    ## Defining the self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        aPar <- min(y) * 0.95
        bPar <- max(y) * 1.05

        ## Linear regression on pseudo y values through origin
        ##  to determine the parameter c
        pseudoY <- log(- ( (y - aPar)/(bPar - aPar) - 1 ) )
        cPar <- coef(lm(pseudoY ~ I(-x) - 1))
        return(c(aPar, bPar, cPar)[notFixed])
    }

    ## Defining names
    names <- names[notFixed]

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        ## Creating the full parameter vector
        ##  containing both estimated and fixed parameters
        parmVec[notFixed] <- parm
        
        ## Converting to relative scale
        if (type == "absolute") 
        {
            p <- 100*((parmVec[2] - respl)/(parmVec[2] - parmVec[1]))
        } else {  
            p <- respl
        }
        
        ## Calculating ED value relative to the control
        if (reference == "control")
        {
            p <- 100 - p
        }
    
        tempVal <- log((100-p)/100)
#        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])
        pProp <- p / 100
        EDp <- -log(pProp)/parmVec[3]

        EDder <- 
#        EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
#        0, 0, 1/parmVec[4], 
#        exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))
        c(0, 0, log(pProp) / (parmVec[3]^2)) 

        return(list(EDp, EDder[notFixed]))
    }
    
    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        log(- ( (y - parmVec[1]) / (parmVec[2] - parmVec[1]) - 1 ) ) / (-parmVec[3])
    }    
    
    ## Returning the function components
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = NULL, deriv2 = NULL, derivx = NULL,
    edfct = edfct, inversion = invfct, 
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName), 
    text = ifelse(missing(fctText), "Asymptotic regression", fctText), 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#' @rdname arandaordaz
#'
#' @title Asymptotic Regression Model
#'
#' @description
#' \code{AR.2} is the two-parameter asymptotic regression model with the lower limit
#' fixed at 0.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names vector of character strings giving the names of the parameters
#'   (should not contain ":").
#'
#' @return A list of class \code{drcMean}.
#'
#' @note The function is for use with \code{\link{drm}}.
#'
#' @examples
#' ## Asymptotic regression on methionine data
#' met.as.m1 <- drm(gain ~ dose, product, data = methionine, fct = AR.3(),
#' pmodels = list(~1, ~factor(product), ~factor(product)))
#' plot(met.as.m1, log = "", ylim = c(1450, 1800))
#' summary(met.as.m1)
#'
#' @keywords models nonlinear
"AR.2" <-
function(fixed = c(NA, NA), names = c("b", "c"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( asymreg(fixed = c(0, fixed[1:2]), 
    names = c("a", names[1:2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Asymptotic regression with lower limit fixed at 0") )
}

#' @rdname arandaordaz
#'
#' @title Asymptotic Regression Model
#'
#' @description
#' \code{AR.3} is the three-parameter asymptotic regression model.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names vector of character strings giving the names of the parameters
#'   (should not contain ":").
#'
#' @return A list of class \code{drcMean}.
#'
#' @keywords models nonlinear
"AR.3" <-
function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( asymreg(fixed, names, 
    fctName = as.character(match.call()[[1]])) )
}
