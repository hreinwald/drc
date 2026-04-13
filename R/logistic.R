#' The general asymmetric five-parameter logistic model
#'
#' The five-parameter logistic model given by the expression
#' \deqn{f(x) = c + \frac{d - c}{(1 + \exp(b(x - e)))^f}}
#'
#' This model differs from the log-logistic in that it uses \code{x} directly
#' rather than \code{log(x)}. It is sometimes referred to as the Boltzmann model.
#'
#' @param fixed numeric vector of length 5. Specifies which parameters are fixed
#'   and at what value they are fixed. \code{NA} indicates that the corresponding
#'   parameter is not fixed.
#' @param names character vector of length 5 giving the names of the parameters
#'   \code{(b, c, d, e, f)}. Default is \code{c("b", "c", "d", "e", "f")}.
#' @param method character string indicating the self starter function to use
#'   (\code{"1"}, \code{"2"}, \code{"3"}, or \code{"4"}).
#' @param ssfct a self starter function to be used. If \code{NULL} (default),
#'   a built-in self starter is selected via \code{method}.
#' @param fctName optional character string used internally to overwrite the
#'   function name.
#' @param fctText optional character string used internally to overwrite the
#'   description text.
#'
#' @return A list of class \code{"Boltzmann"} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{L.3}}, \code{\link{L.4}}, \code{\link{L.5}},
#'   \code{\link{llogistic}}
#'
#' @keywords models nonlinear
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.4())
"logistic" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2]+(parmMat[,3]-parmMat[,2])/((1+exp(parmMat[,1]*(dose-parmMat[,4])))^parmMat[,5])
    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- logistic.ssf(method, fixed)     
    }
   
    ## Defining names
    names <- names[notFixed]

    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(dose - parmMat[, 4]))
        t3 <- (1 + t2)^(2*parmMat[, 5])
        t4 <- parmMat[, 5]*((1 + t2)^(-parmMat[, 5] - 1))
        t5 <- (1 + t2)^(parmMat[, 5])                  

        cbind( -t1*t2*t4*(dose - parmMat[ , 4]), 
               1 - 1/t5, 
               1/t5, 
               t1*t2*t4*parmMat[, 1], 
               -t1*log(1+t2)/t5 )[, notFixed]
    }
        
    deriv2 <- NULL
    
    ##Defining the first derivatives (in x)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        temp1 <- exp(parmMat[, 1]*(x - parmMat[, 4]))       
        
        (-parmMat[, 5]*(parmMat[, 3] - parmMat[, 2])*temp1*parmMat[, 1])/((1 + temp1)^(parmMat[, 5] + 1))
    }

    ## Defining the ED function
    edfct <- function(parm, respl, reference = "control", type = "relative", ...)
    {
        parmVec[notFixed] <- parm

        ## Convert absolute response level to relative.
        ## Note: unlike log-logistic models where b < 0 means decreasing,
        ## the logistic model has b < 0 = increasing.  EDhelper's p-swap
        ## (for b < 0, relative type) would be wrong here, so we perform
        ## only the absolute-to-relative conversion inline.
        if (identical(type, "absolute")) {
            p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
        } else {
            p <- respl
        }
    
        ## deriv(~e + log((100/(100-p))^(1/f) - 1) / b, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
        ## evaluated at the R prompt
        EDderFct <- 
        function (b, c, d, e, f) 
        {
            .expr2 <- 100/p
            .expr4 <- .expr2^(1/f)
            .expr5 <- .expr4 - 1
            .expr6 <- log(.expr5)
            .value <- e + .expr6/b
            .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
            .grad[, "b"] <- -(.expr6/b^2)
            .grad[, "c"] <- 0
            .grad[, "d"] <- 0
            .grad[, "e"] <- 1
            .grad[, "f"] <- -(.expr4 * (log(.expr2) * (1/f^2))/.expr5/b)
            attr(.value, "gradient") <- .grad
            .value
        }
        EDcalc <- EDderFct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDp <- as.numeric(EDcalc)
        EDder <- attr(EDcalc, "gradient")

        ## Fix: correct c and d derivatives for absolute type using central differences.
        ## The analytical derivatives above miss the chain-rule contribution from
        ## the absolute-to-relative conversion, where p depends on c and d.
        if (identical(type, "absolute")) {
            .edval <- function(pv) {
                p0 <- 100 * ((pv[3] - respl) / (pv[3] - pv[2]))
                .expr2 <- 100 / p0
                .expr4 <- .expr2^(1 / pv[5])
                .expr5 <- .expr4 - 1
                pv[4] + log(.expr5) / pv[1]
            }
            .eps <- .Machine$double.eps
            for (.i in c(2, 3)) {
                .h <- if (abs(parmVec[.i]) > sqrt(.eps)) abs(parmVec[.i]) * .eps^(1/3) else .eps^(1/3)
                .pvUp <- replace(parmVec, .i, parmVec[.i] + .h)
                .pvDn <- replace(parmVec, .i, parmVec[.i] - .h)
                EDder[.i] <- (.edval(.pvUp) - .edval(.pvDn)) / (2 * .h)
            }
        }

        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + parmVec[4]
    }

    ## Defining return list
    returnList <- list(fct = fct, ssfct = ssfct, names = names, 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    inversion = invfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Logistic (ED50 as parameter)", fctText),    
    noParm = sum(is.na(fixed)), fixed = fixed)

    class(returnList) <- "Boltzmann"
    invisible(returnList)
}

#' Three-parameter logistic model
#'
#' A three-parameter logistic model with the lower limit fixed at 0, given by
#' \deqn{f(x) = \frac{d}{1 + \exp(b(x - e))}}
#'
#' @param fixed numeric vector of length 3. Specifies which parameters are fixed
#'   and at what value they are fixed. \code{NA} indicates that the corresponding
#'   parameter is not fixed.
#' @param names character vector of length 3 giving the names of the parameters
#'   \code{(b, d, e)}. Default is \code{c("b", "d", "e")}.
#' @param ... additional arguments passed to \code{\link{logistic}}.
#'
#' @return A list of class \code{"Boltzmann"} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{logistic}}, \code{\link{L.4}}, \code{\link{L.5}}
#'
#' @keywords models nonlinear
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.3())
"L.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed[1], 0, fixed[2:3], 1), names = c(names[1], "c", names[2:3], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Logistic (ED50 as parameter) with lower limit fixed at 0", ...))
}

#' Four-parameter logistic model
#'
#' A four-parameter logistic model (symmetric, with \code{f = 1}), given by
#' \deqn{f(x) = c + \frac{d - c}{1 + \exp(b(x - e))}}
#'
#' @param fixed numeric vector of length 4. Specifies which parameters are fixed
#'   and at what value they are fixed. \code{NA} indicates that the corresponding
#'   parameter is not fixed.
#' @param names character vector of length 4 giving the names of the parameters
#'   \code{(b, c, d, e)}. Default is \code{c("b", "c", "d", "e")}.
#' @param ... additional arguments passed to \code{\link{logistic}}.
#'
#' @return A list of class \code{"Boltzmann"} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{logistic}}, \code{\link{L.3}}, \code{\link{L.5}}
#'
#' @keywords models nonlinear
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.4())
"L.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed, 1), names = c(names, "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Logistic (ED50 as parameter)", ...))
}

#' Five-parameter generalized logistic model
#'
#' A five-parameter generalized logistic model (asymmetric when \code{f != 1}),
#' given by
#' \deqn{f(x) = c + \frac{d - c}{(1 + \exp(b(x - e)))^f}}
#'
#' @param fixed numeric vector of length 5. Specifies which parameters are fixed
#'   and at what value they are fixed. \code{NA} indicates that the corresponding
#'   parameter is not fixed.
#' @param names character vector of length 5 giving the names of the parameters
#'   \code{(b, c, d, e, f)}. Default is \code{c("b", "c", "d", "e", "f")}.
#' @param ... additional arguments passed to \code{\link{logistic}}.
#'
#' @return A list of class \code{"Boltzmann"} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{logistic}}, \code{\link{L.3}}, \code{\link{L.4}}
#'
#' @keywords models nonlinear
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.5())
"L.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return(logistic(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), 
    fctText = "Generalised logistic (ED50 as parameter)", ...))
}
