#' The four-parameter Weibull (type 2) model
#'
#' Provides a general framework for the four-parameter Weibull type 2 model
#' given by the equation
#' \deqn{f(x) = c + (d - c)(1 - \exp(-\exp(b(\log(x) - \log(e)))))}
#'
#' @details
#' The \code{method} argument determines how starting values for the parameters
#' \code{b} and \code{e} are estimated (the starting values for \code{c} and
#' \code{d} are always based on the range of the response values). Four methods
#' are available:
#' \describe{
#'   \item{\code{"1"} (default)}{Linear regression on transformed data. Applies a
#'     complementary log-log transformation to the response and a log
#'     transformation to the dose, then fits a linear regression to estimate
#'     starting values for \code{b} and \code{e}.}
#'   \item{\code{"2"}}{Anke's procedure. Estimates \code{e} by finding the dose
#'     at which the response crosses the midpoint between \code{c} and \code{d},
#'     then estimates \code{b} as the median of back-calculated values.}
#'   \item{\code{"3"}}{Stepwise approach. Identifies where the mean response
#'     crosses the midpoint between \code{c} and \code{d} and uses the
#'     corresponding dose as the starting value for \code{e}. The starting value
#'     for \code{b} is based on the sign of the slope at that point.}
#'   \item{\code{"4"}}{Normolle's procedure. Uses the mean of the dose range as
#'     an initial estimate for \code{e}, then estimates \code{b} and \code{e}
#'     using median-based back-calculations.}
#' }
#'
#' @param fixed numeric vector of length 4, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 4 giving the names of the parameters
#'   (default \code{c("b", "c", "d", "e")}).
#' @param method character string indicating the self starter method to use for
#'   obtaining starting values. One of \code{"1"} (default), \code{"2"},
#'   \code{"3"}, or \code{"4"}. See Details.
#' @param ssfct a self starter function. If \code{NULL} (default), a built-in
#'   self starter is used based on \code{method}.
#' @param fctName optional character string used internally for the function name.
#' @param fctText optional character string used internally for the function description.
#'
#' @return A list containing the nonlinear function, self starter function,
#'   and parameter names. The list has class \code{"Weibull-2"}.
#'
#' @author Christian Ritz
#'
#' @references Seber, G. A. F. and Wild, C. J. (1989)
#'   \emph{Nonlinear Regression}, New York: Wiley & Sons (pp. 338--339).
#'
#' @seealso \code{\link{weibull1}}, \code{\link{W2.2}}, \code{\link{W2.3}},
#'   \code{\link{W2.4}}
#'
#' @keywords models nonlinear
"weibull2" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 4
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
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        parmMat[,2] + (parmMat[,3] - parmMat[,2]) * (1 - exp(-exp(parmMat[,1] *(log(dose) - log(parmMat[,4])))))
    }


    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- weibull2.ssf(method, fixed)
    }       

   
    ## Defining names
    w2.names <- names[notFixed]


    ## Defining derivatives
    deriv1 <- function(dose, parm)
              {
                  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
                  parmMat[, notFixed] <- parm

                  t1 <- parmMat[, 3] - parmMat[, 2]
                  t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
                  t3 <- exp(-t2)

                  derMat <- as.matrix(cbind( t1*xexplogx(dose/parmMat[, 4], parmMat[, 1]), 
                                             1 - (1 - t3), 
                                             1 - t3, 
                                             -t1*xexpx(dose/parmMat[, 4], parmMat[, 1])*parmMat[, 1]/parmMat[, 4] ))
                  return(derMat[, notFixed])
              }
    deriv2 <- NULL

 
    ## Defining the first derivative (in x=dose)
    ##  based on deriv(~c+(d-c)*(1 - exp(-exp(b*(log(x)-log(e))))), "x", function(x, b,c,d,e){})
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      .expr1 <- parmMat[, 3] - parmMat[, 2]
      .expr6 <- exp(parmMat[, 1] * (log(x) - log(parmMat[, 4])))
      .expr8 <- exp(-.expr6)
      .value <- parmMat[, 2] + .expr1 * (1 - .expr8)
      .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
      .grad[, "x"] <- .expr1 * (.expr8 * (.expr6 * (parmMat[, 1] * (1/x))))
      .grad
    }
    

    ## Defining the ED function
    edfct <- function(parm, p, reference, type, ...)
    {   
        parmVec[notFixed] <- parm
        respl <- p  # save original response level

        p <- absToRel(parmVec, p, type)

        ## Reversing p
        if (identical(type, "absolute") && (parmVec[1] > 0) && (reference == "control"))
        {
            p <- 100 - p
        }
               
        result <- weibull1(fixed, names)$edfct(parm, p, reference, "relative", ...) 

        ## Fix: correct c and d derivatives for absolute type using central differences.
        ## The delegation to weibull1 with type="relative" produces zero derivatives
        ## for c and d, missing the chain-rule contribution from the
        ## absolute-to-relative conversion (absToRel) where p depends on c and d.
        if (identical(type, "absolute")) {
            .edval <- function(pv) {
                p0 <- absToRel(pv, respl, type)
                # Replicate weibull2's reversal (for b > 0 and absolute type)
                if (pv[1] > 0 && identical(reference, "control")) p0 <- 100 - p0
                # Replicate weibull1's EDhelper swap (for b < 0 and relative type)
                if (pv[1] < 0 && identical(reference, "control")) p0 <- 100 - p0
                tv0 <- log(-log((100 - p0) / 100))
                exp(tv0 / pv[1] + log(pv[4]))
            }
            .eps <- .Machine$double.eps
            .nfIdx <- which(notFixed)
            for (.i in c(2, 3)) {
                if (!notFixed[.i]) next
                .h <- if (abs(parmVec[.i]) > sqrt(.eps)) abs(parmVec[.i]) * .eps^(1/3) else .eps^(1/3)
                .pvUp <- replace(parmVec, .i, parmVec[.i] + .h)
                .pvDn <- replace(parmVec, .i, parmVec[.i] - .h)
                .pos <- which(.nfIdx == .i)
                if (length(.pos) == 1L) {
                    result[[2]][.pos] <- (.edval(.pvUp) - .edval(.pvDn)) / (2 * .h)
                }
            }
        }

        result
    }


    returnList <-
    list(fct = fct, ssfct = ssfct, names = w2.names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName),as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Weibull (type 2)", fctText),     
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "Weibull-2"
    invisible(returnList)
}


#' Two-parameter Weibull (type 2) model
#'
#' A two-parameter Weibull type 2 model with the lower limit fixed at 0 and the
#' upper limit fixed at a specified value. The model is given by the equation
#' \deqn{f(x) = \mathrm{upper} \cdot (1 - \exp(-\exp(b(\log(x) - \log(e)))))}
#' This model is primarily intended for binomial/quantal responses.
#'
#' @param upper numeric value giving the fixed upper limit (default 1).
#' @param fixed numeric vector of length 2, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 2 giving the names of the parameters
#'   (default \code{c("b", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{weibull2}}, \code{\link{W2.3}}, \code{\link{W2.4}},
#'   \code{\link{W1.2}}
#'
#' @examples
#' earthworms.m1 <- drm(number/total ~ dose, weights = total,
#'   data = earthworms, fct = W2.2(), type = "binomial")
#'
#' @keywords models nonlinear
"W2.2" <- function(
upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1], 0, upper, fixed[2]), names = c(names[1], "c", "d", names[2]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Weibull (type 2)", upper), ...))
}

#' Three-parameter Weibull (type 2) model
#'
#' A three-parameter Weibull type 2 model with the lower limit fixed at 0.
#' The model is given by the equation
#' \deqn{f(x) = d \cdot (1 - \exp(-\exp(b(\log(x) - \log(e)))))}
#'
#' @param fixed numeric vector of length 3, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 3 giving the names of the parameters
#'   (default \code{c("b", "d", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{weibull2}}, \code{\link{W2.2}}, \code{\link{W2.4}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2.3())
#'
#' @keywords models nonlinear
"W2.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1], 0, fixed[2:3]), names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 2)"), ...))
}

#' Three-parameter Weibull (type 2) model with upper limit fixed
#'
#' A three-parameter Weibull type 2 model with the upper limit fixed at a
#' specified value. The model is given by the equation
#' \deqn{f(x) = c + (\mathrm{upper} - c)(1 - \exp(-\exp(b(\log(x) - \log(e)))))}
#'
#' @param upper numeric value giving the fixed upper limit (default 1).
#' @param fixed numeric vector of length 3, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 3 giving the names of the parameters
#'   (default \code{c("b", "c", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{weibull2}}, \code{\link{W2.3}}, \code{\link{W2.4}}
#'
#' @keywords models nonlinear
"W2.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1:2], upper, fixed[3]), 
    names = c(names[1:2], "d", names[3]), 
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Weibull (type 2)", upper), ...))
}

#' Four-parameter Weibull (type 2) model
#'
#' A four-parameter Weibull type 2 model. The model is given by the equation
#' \deqn{f(x) = c + (d - c)(1 - \exp(-\exp(b(\log(x) - \log(e)))))}
#'
#' @param fixed numeric vector of length 4, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 4 giving the names of the parameters
#'   (default \code{c("b", "c", "d", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{weibull2}}, \code{\link{W2.2}}, \code{\link{W2.3}}
#'
#' @examples
#' terbuthylazin.m1 <- drm(rgr ~ dose, data = terbuthylazin, fct = W2.4())
#'
#' @keywords models nonlinear
"W2.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}

    return(weibull2(fixed = fixed, names = names,    
    fctName = as.character(match.call()[[1]]),
    fctText = "Weibull (type 2)", ...))
}

#' Two-parameter asymptotic regression model
#'
#' A two-parameter asymptotic regression model where \code{b} is fixed at 1 and
#' the lower limit is fixed at 0. The model is given by the equation
#' \deqn{f(x) = d \cdot (1 - \exp(-x / e))}
#'
#' @param fixed numeric vector of length 2, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 2 giving the names of the parameters
#'   (default \code{c("d", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{AR.3}}, \code{\link{weibull2}}, \code{\link{EXD.2}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.2())
#'
#' @keywords models nonlinear
"AR.2" <-
function(fixed = c(NA, NA), names = c("d", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(1, 0, fixed[1:2]), 
    names = c("b", "c", names[1:2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Asymptotic regression"), ...))
}

#' Three-parameter shifted asymptotic regression model
#'
#' A three-parameter asymptotic regression model where \code{b} is fixed at 1.
#' The model is given by the equation
#' \deqn{f(x) = c + (d - c)(1 - \exp(-x / e))}
#'
#' @param fixed numeric vector of length 3, specifying fixed parameters (use \code{NA} for
#'   parameters that should be estimated).
#' @param names character vector of length 3 giving the names of the parameters
#'   (default \code{c("c", "d", "e")}).
#' @param ... additional arguments passed to \code{\link{weibull2}}, most
#'   notably \code{method} (a character string: \code{"1"} (default),
#'   \code{"2"}, \code{"3"}, or \code{"4"}) which selects the self-starter
#'   method for obtaining starting values. See \code{\link{weibull2}} for
#'   details.
#'
#' @return A list of class \code{"Weibull-2"} as returned by \code{\link{weibull2}}.
#'
#' @seealso \code{\link{AR.2}}, \code{\link{weibull2}}, \code{\link{EXD.3}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.3())
#'
#' @keywords models nonlinear
"AR.3" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(1, fixed[1:3]), 
    names = c("b", names[1:3]), 
    fctName = as.character(match.call()[[1]]),
    fctText = "Shifted asymptotic regression", ...))
}
