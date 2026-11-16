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
#' The parameters are: \code{b}, the rate of the exponential decline above the threshold; \code{c},
#' the lower asymptote (the response approached at high doses); \code{d}, the upper asymptote, i.e.
#' the constant response up to the threshold; and \code{e}, the threshold itself, the
#' no-effect concentration (NEC), which is estimated directly as a model parameter.
#'
#' \strong{NEC versus ED.} The threshold \code{e} (the NEC) and effective doses are different
#' quantities. \code{e} is the dose where the decline starts and is read directly from
#' \code{coef()} (with its standard error from \code{summary()} and interval from
#' \code{confint()}). \code{\link{ED}} returns the dose at which the response has fallen by
#' \eqn{p} percent of the distance between the upper and lower asymptote, i.e. where
#' \eqn{f(x) = d - (p/100)(d - c)}, which gives
#' \deqn{ED_p = e - \log(1 - p/100)/b.}
#' Because the response is flat up to \code{e}, every \eqn{ED_p} with \eqn{0 < p < 100} is larger
#' than the NEC. For example, \code{ED(fit, 10)} is the EC10 and is not the NEC. The ED is
#' available for \code{type = "relative"} and \code{type = "absolute"}; it requires \eqn{b > 0}
#' and, for absolute levels, a response strictly between \code{c} and \code{d}. Otherwise
#' \code{NA} is returned with a warning. Standard errors and intervals from
#' \code{\link{predict}} are based on numerical derivatives.
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
#' ## Fit with drm() and plot; the parameter e is the NEC
#' nec.m1 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
#' summary(nec.m1)
#' plot(nec.m1, type = "all", main = "NEC model")
#'
#' ## NEC (threshold e) versus effective doses: ED() is larger than the NEC
#' coef(nec.m1)["e:(Intercept)"]            # the NEC, estimated directly
#' ED(nec.m1, c(10, 50), display = FALSE)   # EC10 and EC50, both above the NEC
#' ED(nec.m1, 10, interval = "delta", display = FALSE)
#'
#' ## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
#' nec.c0 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, 0, NA, NA)))
#' nec.c3 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, 3, NA, NA)))
#' plot(nec.m1, type = "all", main = "Lower asymptote c")
#' plot(nec.c0, type = "none", add = TRUE, col = "darkgreen")
#' plot(nec.c3, type = "none", add = TRUE, col = "red")
#' legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' ## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
#' nec.d6 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, NA, 6, NA)))
#' nec.d9 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, NA, 9, NA)))
#' plot(nec.m1, type = "all", main = "Upper asymptote d")
#' plot(nec.d6, type = "none", add = TRUE, col = "darkgreen")
#' plot(nec.d9, type = "none", add = TRUE, col = "red")
#' legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
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

    ## Defining the ED function: dose at which the response has dropped by p% of (d - c),
    ## i.e. solves c + (d - c) * exp(-b * (x - e)) = d - p/100 * (d - c) for x > e.
    ## Only defined for b > 0 and 0 < p < 100; otherwise NA.
    edfct <- function(parm, respl, reference = "control", type = "relative", ...)
    {
        parmVec[notFixed] <- parm
        p <- absToRel(parmVec, respl, type)

        if (!is.finite(parmVec[1]) || parmVec[1] <= 0)
        {
            warning("ED is not defined for this NEC fit: the slope parameter b is ",
                    format(parmVec[1], digits = 4), " but must be positive (the response ",
                    "must decay beyond the threshold e). NA returned. Check the fit, e.g. ",
                    "whether the data show a decrease with dose, or refit with different ",
                    "starting values (start =).", call. = FALSE)
            return(list(NA_real_, rep(NA_real_, sum(notFixed))))
        }
        if (!is.finite(p) || p <= 0 || p >= 100)
        {
            if (identical(type, "absolute"))
            {
                msg <- paste0("the absolute response level ", format(respl, digits = 4),
                              " is not strictly between the asymptotes c = ",
                              format(parmVec[2], digits = 4), " and d = ",
                              format(parmVec[3], digits = 4),
                              ". Choose a level inside this range.")
            } else {
                msg <- paste0("the relative response level ", format(respl, digits = 4),
                              " is outside ]0, 100[. Choose a percentage strictly between ",
                              "0 and 100.")
            }
            warning("ED is not defined for NEC models: ", msg, " NA returned.", call. = FALSE)
            return(list(NA_real_, rep(NA_real_, sum(notFixed))))
        }

        logTerm <- log(1 - p/100)
        EDp <- parmVec[4] - logTerm / parmVec[1]
        ## Derivatives wrt b, c, d, e (p treated as constant; for type = "absolute"
        ## ED.drc replaces the gradient by central differences)
        EDder <- c(logTerm / parmVec[1]^2, 0, 0, 1)

        list(EDp, EDder[notFixed])
    }
  
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
#' @examples
#' ## Response relative to the control mean, upper limit fixed at 1
#' ryegrass.rel <- transform(ryegrass, rel = rootl / mean(rootl[conc == 0]))
#' nec.m2 <- drm(rel ~ conc, data = ryegrass.rel, fct = NEC.2(upper = 1))
#' plot(nec.m2, type = "all", main = "NEC.2 (c = 0, d = 1)")
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
#' @examples
#' ## Fit with drm() and plot
#' nec.m3 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.3())
#' plot(nec.m3, type = "all", main = "NEC.3 (lower limit fixed at 0)")
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
#' @examples
#' ## Fit with drm() and plot
#' nec.m4 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
#' plot(nec.m4, type = "all", main = "NEC.4")
#'
#' @keywords models nonlinear
"NEC.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), ...) )
}
