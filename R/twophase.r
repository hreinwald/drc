#' Two-Phase Dose-Response Model
#'
#' A seven-parameter dose-response model combining two log-logistic components,
#' useful for describing more complex dose-response patterns.
#'
#' Following Groot \emph{et al} (1996) the two-phase model function is:
#'
#' \deqn{f(x) = c + \frac{d1-c}{1+\exp(b1(\log(x)-\log(e1)))} + \frac{d2}{1+\exp(b2(\log(x)-\log(e2)))}}
#'
#' For each of the two phases, the parameters have the same interpretation as in
#' the ordinary log-logistic model.
#'
#' @param fixed numeric vector specifying which parameters are fixed and at what value
#'   they are fixed. NAs are used for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters
#'   (should not contain ":"). The default is reasonable.
#' @param fctName optional character string used internally by convenience functions.
#' @param fctText optional character string used internally by convenience functions.
#'
#' @return A list containing the nonlinear function, the self starter function,
#'   and the parameter names.
#'
#' @references Groot, J. C. J., Cone, J. W., Williams, B. A., Debersaques, F. M. A.,
#'   Lantinga, E. A. (1996) Multiphasic analysis of gas production kinetics for
#'   in vitro fermentation of ruminant feeds,
#'   \emph{Animal Feed Science Technology}, \bold{64}, 77--89.
#'
#' @author Christian Ritz
#'
#' @seealso The basic component in the two-phase model is the log-logistic model
#'   \code{\link{llogistic}}.
#'
#' @keywords models nonlinear
"twophase" <- function(
fixed = c(NA, NA, NA, NA, NA, NA, NA), names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2"), 
fctName, fctText)
{   
    ## Checking arguments
    numParm <- 7
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
#        print("A")
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
#        print(notFixed)
#        print(parm)
#        print(parmMat[, notFixed])
        parmMat[, notFixed] <- parm
#        print("B")
#        LL.4(fixed[1:4])$fct(dose, parmMat[, 1:4]) + LL.3(fixed[5:7])$fct(dose, parmMat[, 5:7])

        fixed1.4 <- fixed[1:4]
        fixed5.7 <- fixed[5:7] 
        LL.4()$fct(dose, parmMat[, 1:4, drop = FALSE]) + LL.3()$fct(dose, parmMat[, 5:7, drop = FALSE])
    }

    ## Defining self starter function
    ssfct <- function(dframe)
    {  
#        first4 <- drc:::llogistic.ssf(fixed = fixed[1:4])(dframe)  # drc::: not need
#        first4 <- drc:::llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
        first4 <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
    
#    print(c(first4[1:2], first4[3]/2, first4[4]/3, first4[1], first4[3]/2, first4[4])[is.na(fixed)])

#        c(first4[1:2], first4[3]/2, first4[4]/3, first4[1], first4[3]/2, first4[4])[is.na(fixed)]
        c(first4[1:2], first4[3]/2, first4[4]/3, first4[1], first4[3], first4[4])[is.na(fixed)]
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
    text = ifelse(missing(fctText), "Two-phase", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "two-phase"
    invisible(returnList)
}
