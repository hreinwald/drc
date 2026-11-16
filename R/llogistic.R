#' The log-logistic function
#'
#' A very general way of specifying log-logistic models under various
#' constraints on parameters.
#'
#' The five-parameter log-logistic function is given by the expression
#' \deqn{f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
#'
#' @details
#' The parameters are: \code{b}, the slope around \code{e} (positive for decreasing curves, negative
#' for increasing curves); \code{c}, the lower asymptote; \code{d}, the upper asymptote (the response
#' of the untreated control for a decreasing curve); \code{e}, the ED50, i.e. the dose half-way between
#' \code{c} and \code{d} when \code{f = 1}; and \code{f}, the asymmetry parameter (\code{f = 1} gives the
#' symmetric model). Fixing parameters via \code{fixed} yields the convenience models
#' \code{\link{LL.2}}, \code{\link{LL.3}}, \code{\link{LL.3u}}, \code{\link{LL.4}} and \code{\link{LL.5}}.
#'
#' @param fixed numeric vector of length 5, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 5, specifying the names of the
#'   parameters: b, c, d, e, f.
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used.
#' @param fctName optional character string used internally.
#' @param fctText optional character string used internally.
#'
#' @return A list containing the nonlinear function, the self starter function,
#'   and the parameter names.
#'
#' @author Christian Ritz
#'
#' @references
#'   Finney, D. J. (1979).
#'
#'   Seber, G. A. F. and Wild, C. J. (1989).
#'
#' @seealso \code{\link{LL.2}}, \code{\link{LL.3}}, \code{\link{LL.4}},
#'   \code{\link{LL.5}}
#'
#' @examples
#' ## Four-parameter log-logistic model (f fixed at 1): fit with drm() and plot
#' m1 <- drm(rootl ~ conc, data = ryegrass, fct = llogistic(fixed = c(NA, NA, NA, NA, 1)))
#' plot(m1, type = "all", main = "Four-parameter log-logistic model")
#'
#' ## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
#' m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 0, NA, NA)))
#' m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 3, NA, NA)))
#' plot(m1, type = "all", main = "Lower asymptote c")
#' plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
#' plot(m.c3, type = "none", add = TRUE, col = "red")
#' legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' ## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
#' m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 6, NA)))
#' m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 9, NA)))
#' plot(m1, type = "all", main = "Upper asymptote d")
#' plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
#' plot(m.d9, type = "none", add = TRUE, col = "red")
#' legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' @keywords models nonlinear
"llogistic" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL, 
fctName, fctText)
{
    ## Matching 'adjust' argument
    method <- match.arg(method)
    
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the model function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose/parmMat[, 4]))))^parmMat[, 5])
    }
       
    ## Defining the model function adjusted for scaling
    retFct <- function(doseScaling, respScaling)
    {          
        fct <- function(dose, parm) 
        {
            parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), 
                              nrow(parm), numParm, byrow = TRUE)
            parmMat[, notFixed] <- parm
        
            cParm <- parmMat[, 2]
            cParm + (parmMat[, 3] - cParm)/((1 + exp(parmMat[, 1]*(log(dose / parmMat[, 4]))))^parmMat[, 5])
        }
        fct        
    }

    ## Defining the derivative in x adjusted for scaling
    retFctDx <- function(doseScaling, respScaling)
    {          
      fct <- function(dose, parm) 
      {
        parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), 
                          nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        temp1 <- dose/parmMat[, 4]
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec <- (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec
      }
      fct        
    }
    
        
    ## Defining the scale function
    scaleFct <- function(doseScaling, respScaling)
    {        
        c(1, respScaling, respScaling, doseScaling, 1)[notFixed]
    }    

    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- llogistic.ssf(method, fixed)
    }
   
    ## Defining names
    names <- names[notFixed]
    
    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t5 <- (1 + t2)^parmMat[, 5]                  

        cbind( -t1 * xlogx(dose/parmMat[, 4], parmMat[, 1], parmMat[, 5] + 1) * parmMat[, 5],
               1 - 1/t5, 
               1/t5, 
               t1 * parmMat[, 5] * divAtInf(t2, (1 + t2)^(parmMat[, 5] + 1)) * parmMat[, 1] / parmMat[, 4], 
               -t1 * divAtInf(log(1+t2), t5) )[, notFixed]
    }
        
    deriv2 <- NULL


    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        temp1 <- x/parmMat[, 4]
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec <- (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec
    }


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)

        tempVal <- log((100-p)/100)
        expTerm <- exp(-tempVal/parmVec[5])

        # Check if expTerm - 1 is valid (must be positive for log)
        # Handle NaN from tempVal or expTerm being invalid
        if (is.na(expTerm) || expTerm <= 1) {
            # ED value is outside the valid range or model is ill-conditioned
            EDp <- Inf
            EDder <- rep(NA, 5)
        } else {
            EDp <- parmVec[4]*(expTerm-1)^(1/parmVec[1])

            EDder <-
            EDp*c(-log(expTerm-1)/(parmVec[1]^2),
            0, 0, 1/parmVec[4],
            expTerm*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((expTerm-1)^(-1)))

            ## Fix: correct c and d derivatives for absolute type using central differences.
            ## The analytical derivatives above miss the chain-rule contribution from
            ## the absolute-to-relative conversion (EDhelper), where p depends on c and d.
            if (identical(type, "absolute")) {
                .edval <- function(pv) {
                    p0 <- EDhelper(pv, respl, reference, type)
                    tv0 <- log((100 - p0) / 100)
                    et0 <- exp(-tv0 / pv[5])
                    if (is.na(et0) || et0 <= 1) return(Inf)
                    pv[4] * (et0 - 1)^(1 / pv[1])
                }
                .eps <- .Machine$double.eps
                for (.i in c(2, 3)) {
                    .h <- if (abs(parmVec[.i]) > sqrt(.eps)) abs(parmVec[.i]) * .eps^(1/3) else .eps^(1/3)
                    .pvUp <- replace(parmVec, .i, parmVec[.i] + .h)
                    .pvDn <- replace(parmVec, .i, parmVec[.i] - .h)
                    EDder[.i] <- (.edval(.pvUp) - .edval(.pvDn)) / (2 * .h)
                }
            }
        }

        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm) 
    {
        parmVec[notFixed] <- parm
        
        exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))
    } 
    
    ## Defining functions returning lower and upper limit and monotonicity
    lowerAs <- pickParm(parmVec, notFixed, 2)
    upperAs <- pickParm(parmVec, notFixed, 3)
    monoton <- monoParm(parmVec, notFixed, 1, -1)

    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, inversion = invfct, scaleFct = scaleFct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Log-logistic (ED50 as parameter)", fctText), 
    noParm = sum(is.na(fixed)), lowerAs = lowerAs, upperAs = upperAs, monoton = monoton,
    retFct = retFct, fixed = fixed, retFctDx = retFctDx)
    
    class(returnList) <- "llogistic"
    invisible(returnList)
}

#' Two-parameter log-logistic function
#'
#' A two-parameter log-logistic function with lower limit fixed at 0 and
#' upper limit fixed (default 1), primarily for use with binomial/quantal
#' dose-response data.
#'
#' The two-parameter log-logistic function is given by the expression
#' \deqn{f(x) = \frac{upper}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' @param upper numeric value, the fixed upper limit (default 1).
#' @param fixed numeric vector of length 2, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 2, specifying the names of the
#'   parameters (default: b, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{LL.3}}, \code{\link{LL.4}}, \code{\link{LL.5}},
#'   \code{\link{llogistic}}
#'
#' @examples
#' ## Binomial (quantal) response: fit with drm() and plot
#' earthworms.m1 <- drm(number/total ~ dose, weights = total,
#'   data = earthworms, fct = LL.2(), type = "binomial")
#' plot(earthworms.m1, main = "LL.2 (upper limit fixed at 1)")
#'
#' @keywords models nonlinear
"LL.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, upper, fixed[2], 1), 
    names = c(names[1], "c", "d", names[2], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Log-logistic (ED50 as parameter)", upper),
    ...) )
}

#' @rdname LL.2
l2 <- LL.2

#' Three-parameter log-logistic function
#'
#' A three-parameter log-logistic function with lower limit fixed at 0.
#'
#' The three-parameter log-logistic function is given by the expression
#' \deqn{f(x) = \frac{d}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' @param fixed numeric vector of length 3, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 3, specifying the names of the
#'   parameters (default: b, d, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{LL.2}}, \code{\link{LL.4}}, \code{\link{LL.5}},
#'   \code{\link{llogistic}}
#'
#' @examples
#' ## Fit with drm() and plot
#' ryegrass.model1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
#' plot(ryegrass.model1, type = "all", main = "LL.3 (lower limit fixed at 0)")
#'
#' @keywords models nonlinear
"LL.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1], 0, fixed[2:3], 1), 
    names = c(names[1], "c", names[2:3], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Log-logistic (ED50 as parameter)"), 
    ...) )
}

#' @rdname LL.3
l3 <- LL.3

#' Three-parameter log-logistic function with upper limit fixed
#'
#' A three-parameter log-logistic function with upper limit fixed (default 1),
#' primarily for use with binomial/quantal dose-response data.
#'
#' The three-parameter log-logistic function with upper limit fixed is given by
#' \deqn{f(x) = c + \frac{upper-c}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' @param upper numeric value, the fixed upper limit (default 1).
#' @param fixed numeric vector of length 3, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 3, specifying the names of the
#'   parameters (default: b, c, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{LL.2}}, \code{\link{LL.3}}, \code{\link{LL.4}},
#'   \code{\link{llogistic}}
#'
#' @examples
#' ## Response relative to the control mean, upper limit fixed at 1
#' ryegrass.rel <- transform(ryegrass, rel = rootl / mean(rootl[conc == 0]))
#' ryegrass.m1 <- drm(rel ~ conc, data = ryegrass.rel, fct = LL.3u(upper = 1))
#' plot(ryegrass.m1, type = "all", main = "LL.3u (upper limit fixed at 1)")
#'
#' @keywords models nonlinear
"LL.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed[1:2], upper, fixed[3], 1), 
    names = c(names[1:2], "d", names[3], "f"),
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Log-logistic (ED50 as parameter)", upper), 
    ...) )
}

#' @rdname LL.3u
l3u <- LL.3u

#' Four-parameter log-logistic function
#'
#' A four-parameter log-logistic function.
#'
#' The four-parameter log-logistic function is given by the expression
#' \deqn{f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' @param fixed numeric vector of length 4, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 4, specifying the names of the
#'   parameters (default: b, c, d, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz and Jens C. Streibig
#'
#' @seealso \code{\link{LL.3}}, \code{\link{LL.5}}, \code{\link{llogistic}}
#'
#' @examples
#' ## Several curves fitted with drm() and plotted
#' spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
#' plot(spinach.m1, main = "LL.4, one curve per CURVE level")
#'
#' ## Effect of the asymptotes on a single curve (ryegrass): fix c at 0 / d at 6
#' m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 0, NA, NA)))
#' m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 6, NA)))
#' plot(m1, type = "all", main = "LL.4: free and fixed asymptotes")
#' plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
#' plot(m.d6, type = "none", add = TRUE, col = "red")
#' legend("topright", c("free", "c = 0", "d = 6"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' @keywords models nonlinear
"LL.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(fixed, 1), names = c(names, "f"),
    fctName = as.character(match.call()[[1]]), ...) )
}

#' @rdname LL.4
l4 <- LL.4

#' Five-parameter log-logistic function
#'
#' A five-parameter (generalized) log-logistic function. The function is
#' asymmetric when f differs from 1.
#'
#' The five-parameter log-logistic function is given by the expression
#' \deqn{f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
#'
#' @param fixed numeric vector of length 5, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 5, specifying the names of the
#'   parameters (default: b, c, d, e, f).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{LL.3}}, \code{\link{LL.4}}, \code{\link{llogistic}}
#'
#' @examples
#' ## Fit with drm() and plot; f allows for an asymmetric curve
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.5())
#' plot(ryegrass.m1, type = "all", main = "LL.5")
#'
#' ## Effect of the lower asymptote c
#' ryegrass.c3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.5(fixed = c(NA, 3, NA, NA, NA)))
#' plot(ryegrass.c3, type = "none", add = TRUE, col = "red")
#'
#' @keywords models nonlinear
"LL.5" <-
function(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    return( llogistic(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]),
    fctText = "Generalized log-logistic (ED50 as parameter)", ...) )
}

#' @rdname LL.5
l5 <- LL.5

#' Two-parameter Michaelis-Menten function
#'
#' A two-parameter Michaelis-Menten function where b is fixed at -1, c at 0,
#' and f at 1. Commonly used for enzyme kinetics and weed density studies.
#'
#' The two-parameter Michaelis-Menten function is
#' \deqn{f(x) = \frac{d \cdot x}{e + x}}
#' which is equivalent to \eqn{d/(1+(e/x))}.
#'
#' @param fixed numeric vector of length 2, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 2, specifying the names of the
#'   parameters (default: d, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{MM.3}}, \code{\link{AR.2}}, \code{\link{AR.3}}
#'
#' @examples
#' ## Fit with drm() and plot
#' met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.2())
#' plot(met.mm.m1, main = "MM.2 (Michaelis-Menten)")
#'
#' @keywords models nonlinear
"MM.2" <-
function(fixed = c(NA, NA), names = c("d", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(-1, 0, fixed[1:2], 1), names = c("b", "c", names[1:2], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Michaelis-Menten", 
    ...) )
}

#' Three-parameter Michaelis-Menten function
#'
#' A three-parameter (shifted) Michaelis-Menten function where b is fixed
#' at -1 and f at 1.
#'
#' The three-parameter Michaelis-Menten function is
#' \deqn{f(x) = c + \frac{d-c}{1+(e/x)}}
#'
#' @param fixed numeric vector of length 3, specifying fixed parameters
#'   (use NA for non-fixed parameters).
#' @param names character vector of length 3, specifying the names of the
#'   parameters (default: c, d, e).
#' @param ... additional arguments to \code{\link{llogistic}}.
#'
#' @return See \code{\link{llogistic}}.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{MM.2}}, \code{\link{AR.2}}, \code{\link{AR.3}}
#'
#' @examples
#' ## Fit with drm() and plot
#' met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.3())
#' plot(met.mm.m1, main = "MM.3 (shifted Michaelis-Menten)")
#'
#' @keywords models nonlinear
"MM.3" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( llogistic(fixed = c(-1, fixed[1:3], 1), names = c("b", names[1:3], "f"),
    fctName = as.character(match.call()[[1]]), 
    fctText = "Shifted Michaelis-Menten", 
    ...) )
}
