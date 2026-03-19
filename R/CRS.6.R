#' Generalised Cedergreen-Ritz-Streibig Model for Hormesis
#'
#' A six-parameter extension of the Cedergreen-Ritz-Streibig model for
#' describing hormesis, where the alpha parameter is estimated rather than fixed.
#'
#' The model function is:
#'
#' \deqn{f(x) = c + \frac{d-c+f \exp(-1/x^g)}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' This generalises the five-parameter \code{\link{CRS.5a}} model by estimating
#' the alpha exponent (parameter \eqn{g}) instead of fixing it.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters
#'   (should not contain ":").
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used (optional).
#'
#' @return A list containing the nonlinear model function, the self starter function,
#'   and the parameter names.
#'
#' @author Christian Ritz
#'
#' @note This function is for use with \code{\link{drm}}.
#'
#' @seealso \code{\link{CRS.5a}}, \code{\link{cedergreen}}
#'
#' @keywords models nonlinear
"CRS.6" <- function(
fixed = c(NA, NA, NA, NA, NA, NA), 
names = c("b", "c", "d", "e", "f", "g"),
method = c("1", "2", "3", "4"), 
ssfct = NULL ){
    ## Checking arguments
    numParm <- 6
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2] + (parmMat[,3] - parmMat[,2] + parmMat[,5]*exp(-1/(dose^parmMat[,6])))/(1 + exp(parmMat[,1]*(log(dose) - log(parmMat[,4]))))
    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {   
        ssfct <- function(dframe)
        {
            initval <- c(llogistic()$ssfct(dframe), 0)   
            initval[5] <- (2*(median(dframe[, 2])-initval[2])-(initval[3]-initval[2]))*exp(1/(initval[4]^initval[6]))
    
            return(initval[notFixed])
        }        
    }
 
    ## Defining names
    names <- names[notFixed]

    ## Specifying the derivatives    
    deriv1 <- NULL
    deriv2 <- NULL


    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function    
    edfct <- NULL
    
    ## Finding the maximal hormesis
    maxfct <- NULL
    
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, # lowerc=lowerLimits, upperc=upperLimits, 
    edfct = edfct, maxfct = maxfct,
    name = "CRS.6",
    text = "Generalised Cedergreen-Ritz-Streibig (hormesis)", 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "cedergreen.extended"
    invisible(returnList)
}

