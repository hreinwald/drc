#' The Baroreflex Five-Parameter Dose-Response Model
#'
#' \code{baro5} provides the five-parameter baroreflex model function, allowing
#' specification under various parameter constraints. The model accommodates
#' asymmetric dose-response curves.
#'
#' The five-parameter function is given by:
#'
#' \deqn{y = c + \frac{d-c}{1+f\exp(b1(\log(x)-\log(e))) + (1-f)\exp(b2(\log(x)-\log(e)))}}
#'
#' \deqn{f = 1/(1 + \exp((2b1 b2/|b1+b2|)(\log(x)-\log(e))))}
#'
#' If the difference between b1 and b2 is nonzero, the function is asymmetric.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters
#'   (should not contain ":"). The order is: b1, b2, c, d, e.
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used.
#'
#' @return A list containing the nonlinear model function, the self starter function,
#'   and the parameter names.
#'
#' @references Ricketts, J. H. and Head, G. A. (1999)
#'   A five-parameter logistic equation for investigating asymmetry of curvature
#'   in baroreflex studies.
#'   \emph{Am. J. Physiol. (Regulatory Integrative Comp. Physiol. 46)}, \bold{277}, 441--454.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"baro5" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b1", "b2", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
    
        c <- 2*parmMat[, 1]*parmMat[, 2]/abs(parmMat[, 1]+parmMat[, 2])

        tempVal <- log(dose) - log(parmMat[, 5])
        f <- 1/(1+exp(c*tempVal))
        g <- exp(parmMat[, 1]*tempVal)
        h <- exp(parmMat[, 2]*tempVal)
        parmMat[, 3]+((parmMat[,4]-parmMat[,3])/(1+f*g+(1-f)*h))

    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {   
        ssfct <- function(dframe)
        {
            initval <- (llogistic()$ssfct(dframe))[c(1, 1:4)]   
    
            return(initval[notFixed])
        }        
    }
   
    ## Defining names
    names <- names[notFixed]

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL    
    
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct=edfct, sifct=sifct,
    name = "baro5",
    text = "Baroreflex", 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "baro5"
    invisible(returnList)
}
