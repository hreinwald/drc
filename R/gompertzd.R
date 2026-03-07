#' @title Derivative of the Gompertz function
#'
#' @description
#' \code{gompertzd} provides a way of specifying the derivative of the Gompertz function
#' as a dose-response model.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed.
#'   NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters (should not contain ":").
#'   The default is (notice the order): a, b.
#'
#' @details
#' The derivative of the Gompertz function is defined as
#' \deqn{f(x) = a \exp(bx-a/b(\exp(bx)-1))}
#' For \eqn{a>0} and \eqn{b} not 0, the function is decreasing, equaling \eqn{a} at \eqn{x=0}
#' and approaching 0 at plus infinity.
#'
#' @return A list containing the model function, the self starter function
#'   and the parameter names.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{gompertz}}, \code{\link{drm}}
#'
#' @keywords models nonlinear
"gompertzd" <- function(
fixed = c(NA, NA), names = c("a", "b"))
{   
    ## Checking arguments
    numParm <- 2
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
        
        innerT1 <- parmMat[, 2]*dose
        innerT2 <- parmMat[, 1]/parmMat[, 2]*(exp(innerT1) - 1)
        parmMat[, 1]*exp(innerT1 - innerT2)
    }

    ## Defining the self starter function
    ssfct <- function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]    
    
        aVal <- max(y)
        bVal <- 1
        
        return(c(aVal, bVal))
    }
   
    ## Defining names
    names <- names[notFixed]
        
    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        help1 <- fct(dose, parm)
        help2 <- exp(parmMat[, 2]*dose)
        help3 <- help2 - 1
        help4 <- parmMat[, 1]/parmMat[, 2]
        
        deriva <- help1*(1/parmMat[, 1] - help3/parmMat[, 2])
        derivb <- help1*(dose + help4*help3/parmMat[, 2] + help4*help2*dose)
        
        cbind(deriva, derivb)[, notFixed]       
    }
        
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        fct(x, parm)*(parmMat[, 2] - parmMat[, 1]*exp(parmMat[, 2]*x))
    }

    ## Defining the ED function
    edfct <- NULL
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct, 
    name = "gompertzd",
    text = "Gompertz derivative", 
    noParm = sum(is.na(fixed)))
    
    class(returnList) <- "gompertzd"
    invisible(returnList)
}
