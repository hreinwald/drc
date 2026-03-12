#' @title No Effect Concentration (NEC) dose-response model
#'
#' @description
#' The NEC model is a dose-response model with a threshold below which the response is assumed
#' constant and equal to the control response. It has been proposed as an alternative to both the
#' classical NOEC and the regression-based EC/ED approach.
#'
#' @param fixed numeric vector specifying which parameters are fixed and at what value they are fixed.
#'   NAs are used for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters (should not contain ":").
#'   The default is reasonable (see under 'Usage').
#' @param fctName optional character string used internally by convenience functions.
#' @param fctText optional character string used internally by convenience functions.
#'
#' @details
#' The NEC model function proposed by Pires et al (2002) is:
#' \deqn{f(x) = c + (d-c)\exp(-b(x-e)I(x-e))}
#' where \eqn{I(x-e)} is the indicator function equal to 0 for \eqn{x<=e} and 1 for \eqn{x>e}.
#'
#' @return A list containing the nonlinear function, the self starter function
#'   and the parameter names.
#'
#' @references
#' Pires, A. M., Branco, J. A., Picado, A., Mendonca, E. (2002)
#' Models for the estimation of a 'no effect concentration',
#' \emph{Environmetrics}, \bold{13}, 15--27.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{NEC.2}}, \code{\link{NEC.3}}, \code{\link{NEC.4}}, \code{\link{drm}}
#'
#' @examples
#' nec.m1 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
#' summary(nec.m1)
#' plot(nec.m1)
#'
#' @keywords models nonlinear
"NEC" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
fctName, fctText)
{   
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        cParm <- parmMat[, 2]
        doseDiff <- dose - parmMat[, 4]
        cParm + (parmMat[, 3] - cParm) * exp(-parmMat[, 1] * doseDiff * (doseDiff > 0) )
    }

    ## Defining self starter function
    ssfct <- function(dframe)
    {  
        LLinit <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)  # drc::: not needed      
        c(LLinit[1:3], LLinit[3]/3)[is.na(fixed)]
    }     
      
    ##Defining the first and second derivative (in the parameters) 
    deriv1 <- NULL
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- NULL

    ## Defining the ED function
    edfct <- NULL
  
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "NEC", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "NEC"
    invisible(returnList)
}


#' @title Two-parameter NEC model
#'
#' @description
#' Convenience function for the NEC model with lower limit fixed at 0 and upper limit fixed.
#'
#' @param upper numeric value. The fixed upper limit in the model. Default is 1.
#' @param fixed numeric vector of length 2 specifying fixed parameters (NAs for free parameters).
#' @param names character vector of parameter names.
#' @param ... additional arguments passed to \code{\link{NEC}}.
#'
#' @return A list (see \code{\link{NEC}}).
#'
#' @seealso \code{\link{NEC}}, \code{\link{NEC.3}}, \code{\link{NEC.4}}
#'
#' @keywords models nonlinear
"NEC.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = c(fixed[1], 0, upper, fixed[2]), 
    names = c(names[1], "c", "d", names[2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("NEC", upper),
    ...) )
}

#' @title Three-parameter NEC model
#'
#' @description
#' Convenience function for the NEC model with the lower limit fixed at 0.
#'
#' @param fixed numeric vector of length 3 specifying fixed parameters (NAs for free parameters).
#' @param names character vector of parameter names.
#' @param ... additional arguments passed to \code{\link{NEC}}.
#'
#' @return A list (see \code{\link{NEC}}).
#'
#' @seealso \code{\link{NEC}}, \code{\link{NEC.2}}, \code{\link{NEC.4}}
#'
#' @keywords models nonlinear
"NEC.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = c(fixed[1], 0, fixed[2:3]), 
    names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("NEC"), 
    ...) )
}

#' @title Four-parameter NEC model
#'
#' @description
#' Convenience function for the full four-parameter NEC model.
#'
#' @param fixed numeric vector of length 4 specifying fixed parameters (NAs for free parameters).
#' @param names character vector of parameter names.
#' @param ... additional arguments passed to \code{\link{NEC}}.
#'
#' @return A list (see \code{\link{NEC}}).
#'
#' @seealso \code{\link{NEC}}, \code{\link{NEC.2}}, \code{\link{NEC.3}}
#'
#' @keywords models nonlinear
"NEC.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), ...) )
}
