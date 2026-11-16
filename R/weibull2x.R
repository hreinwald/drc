#' Five-parameter Weibull type 2 model with lag time
#'
#' A five-parameter Weibull type 2 model extended with a lag time parameter
#' \code{t0}. The model is given by the expression
#' \deqn{f(x) = c + (d - c)(1 - \exp(-\exp(b(\log(x - t0) - \log(e)))))}
#' for \eqn{x > t0} and \eqn{f(x) = c} otherwise.
#'
#' The lag time parameter \code{t0} cannot be fixed.
#'
#' @details
#' The parameters are those of \code{\link{weibull2}} (\code{b} slope, \code{c} lower asymptote,
#' \code{d} upper asymptote, \code{e} scale) plus the lag time \code{t0}: the response stays at the
#' lower asymptote \code{c} until the dose exceeds \code{t0}. Effective doses (\code{\link{ED}}) are
#' shifted by \code{t0}; the derivative with respect to \code{t0} (equal to 1) is included in the
#' standard errors.
#'
#' @param fixed numeric vector of length 5. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated (default is \code{rep(NA, 5)}).
#' @param names character vector of length 5 giving the names of the
#'   parameters (default is \code{c("b", "c", "d", "e", "t0")}).
#' @param method character string indicating the self starter method to use.
#'   One of \code{"1"}, \code{"2"}, \code{"3"}, or \code{"4"}.
#' @param ssfct a self starter function. If \code{NULL} (default), a built-in
#'   self starter is used.
#' @param fctName optional character string specifying the function name
#'   (used internally).
#' @param fctText optional character string specifying the function description
#'   (used internally).
#'
#' @return A list of class \code{"Weibull-2"} containing the nonlinear
#'   function, self starter function, and parameter names.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{weibull2}}, \code{\link{W2x.3}}, \code{\link{W2x.4}}
#'
#' @examples
#' ## Fit with drm() and plot (b fixed at 1)
#' m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2x.4())
#' plot(m1, type = "all", main = "Weibull type 2 model with lag time")
#'
#' ## Effect of the lower asymptote c: c free or fixed at 3
#' m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = W2x.4(fixed = c(3, NA, NA, NA)))
#' plot(m1, type = "all", main = "Lower asymptote c")
#' plot(m.c3, type = "none", add = TRUE, col = "red")
#' legend("topright", c("c free", "c = 3"), col = c("black", "red"), lty = 1, bty = "n")
#'
#' @keywords models nonlinear
"weibull2x" <- function(
fixed = rep(NA, 5), names = c("b", "c", "d", "e", "t0"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    
    
    if (!is.na(fixed[5])) {stop("The lag time cannot be fixed")}
    
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        ifelse(dose > parmMat[,5],parmMat[,2] + (parmMat[,3] - parmMat[,2]) * (1 - exp(-exp(parmMat[,1] *(log(dose - parmMat[,5]) - log(parmMat[,4]))))), parmMat[,2]) 
    }


    ## Defining the self starter function  
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- function(dataFrame)
        {
            c((weibull2.ssf(method, head(fixed, -1)))(dataFrame), fivenum(dataFrame[, 1])[2])
        }
    }       

   
    ## Defining names
    w2.names <- names[notFixed]

    
    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL
    derivx <- NULL


    ## Defining the ED function
    edfct <- function(parm, p, reference, type, ...)
    {   
        parmVec[notFixed] <- parm
        
        if (identical(type, "absolute"))
        {
            ## same absolute -> relative conversion as the lag-free Weibull type 2 model
            edResult <- weibull2(head(fixed, -1), head(names, -1))$edfct(head(parm, -1), p, reference, type, ...)
        } else {
            if ( (parmVec[1] > 0) && (reference == "control") ) 
            {
                p <- 100 - p
                reference <- "upper"  # to avoid resetting of p in 'weibull1'    
            }
            if ( (parmVec[1] < 0) && (reference == "control") ) 
            {
                p <- 100 - p
            }
                    
            edResult <- weibull1(head(fixed, -1), head(names, -1))$edfct(head(parm, -1), 100 - p, reference, type, ...)
        }
        # t0 is estimated: its derivative (1) is appended below
        edResult[[1]] <- edResult[[1]] + parmVec[5]
        edResult[[2]] <- c(edResult[[2]], 1)  # d ED / d t0 = 1 (t0 is never fixed)
        edResult
    }
    ## using head( , -1) to remove the lag-time parameter


    ## Returning list of functions and values
    returnList <-
    list(fct = fct, ssfct = ssfct, names = w2.names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName),as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Weibull (type 2)", fctText),     
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "Weibull-2"
    invisible(returnList)
}


#' Three-parameter Weibull type 2 model with lag time
#'
#' A three-parameter Weibull type 2 model with lag time, where \code{b} is
#' fixed at 1 and \code{c} is fixed at 0. This is a convenience wrapper
#' around \code{\link{weibull2x}}.
#'
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated (default is \code{c(NA, NA, NA)}).
#' @param names character vector of length 3 giving the names of the
#'   parameters (default is \code{c("d", "e", "t0")}).
#' @param \dots additional arguments passed to \code{\link{weibull2x}}.
#'
#' @return A list of class \code{"Weibull-2"} containing the nonlinear
#'   function, self starter function, and parameter names.
#'
#' @seealso \code{\link{weibull2x}}, \code{\link{W2x.4}}, \code{\link{W2.3}}
#'
#' @examples
#' ## Fit with drm() and plot
#' spinach.m1 <- drm(SLOPE ~ DOSE, data = spinach, fct = W2x.3())
#' summary(spinach.m1)
#' plot(spinach.m1, main = "W2x.3")
#'
#' @keywords models nonlinear
"W2x.3" <-
function(fixed = c(NA, NA, NA), names = c("d", "e", "t0"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2x(fixed = c(1, 0, fixed[1:3]), names = c("b", "c", names[1:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 2)"), ...))
}            

#' Four-parameter Weibull type 2 model with lag time
#'
#' A four-parameter Weibull type 2 model with lag time, where \code{b} is
#' fixed at 1. This is a convenience wrapper around \code{\link{weibull2x}}.
#'
#' @param fixed numeric vector of length 4. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that should be
#'   estimated (default is \code{c(NA, NA, NA, NA)}).
#' @param names character vector of length 4 giving the names of the
#'   parameters (default is \code{c("c", "d", "e", "t0")}).
#' @param \dots additional arguments passed to \code{\link{weibull2x}}.
#'
#' @return A list of class \code{"Weibull-2"} containing the nonlinear
#'   function, self starter function, and parameter names.
#'
#' @seealso \code{\link{weibull2x}}, \code{\link{W2x.3}}, \code{\link{W2.4}}
#'
#' @examples
#' ## Fit with drm() and plot
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2x.4())
#' plot(ryegrass.m1, type = "all", main = "W2x.4")
#' ED(ryegrass.m1, 50)
#'
#' @keywords models nonlinear
"W2x.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("c", "d", "e", "t0"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2x(fixed = c(1, fixed[1:4]), names = c("b", names[1:4]),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Weibull (type 2)", ...))
}
