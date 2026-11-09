#' Gamma Dose-Response Model
#'
#' A four-parameter dose-response model derived from the cumulative distribution
#' function of the gamma distribution. Only suitable for increasing dose-response data.
#'
#' Following Wheeler and Bailer (2009) the model function is:
#'
#' \deqn{f(x) = c + (d-c) \cdot \mathrm{pgamma}(b \cdot x, e, 1)}
#'
#' The parameters are: \code{b}, the rate (scale) of the dose; \code{c}, the lower asymptote (the
#' response at dose 0); \code{d}, the upper asymptote (approached at high doses); and \code{e}, the
#' shape of the gamma distribution. Effective doses are not available for this model;
#' \code{\link{predict}} uses the analytic derivatives for standard errors.
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
#' @references Wheeler, M. W., Bailer, A. J. (2009)
#'   Comparing model averaging with other model selection strategies for benchmark
#'   dose estimation, \emph{Environmental and Ecological Statistics}, \bold{16}, 37--51.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{drm}}, \code{\link{multi2}}
#'
#' @examples
#' ## Simulated increasing dose-response data: fit with drm() and plot
#' dose <- rep(c(0.5, 1, 2, 4, 8, 16, 32), each = 3)
#' mu <- gammadr()$fct(dose, matrix(c(0.2, 1, 9, 2.5), length(dose), 4, byrow = TRUE))
#' dat <- data.frame(dose = dose, response = mu + 0.25 * sin(seq_along(dose) * 2.3))
#' m1 <- drm(response ~ dose, data = dat, fct = gammadr())
#' plot(m1, type = "all", main = "Gamma dose-response model")
#'
#' ## Effect of the lower asymptote c and the upper asymptote d
#' m.c0 <- drm(response ~ dose, data = dat, fct = gammadr(fixed = c(NA, 0, NA, NA)))
#' m.d7 <- drm(response ~ dose, data = dat, fct = gammadr(fixed = c(NA, NA, 7, NA)))
#' plot(m1, type = "all", main = "Fixed asymptotes")
#' plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
#' plot(m.d7, type = "none", add = TRUE, col = "red")
#' legend("bottomright", c("free", "c = 0", "d = 7"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' @keywords models nonlinear
"gammadr" <- function(
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
        cParm + (parmMat[, 3] - cParm) * pgamma(parmMat[, 1] * dose, parmMat[, 4], 1)
    }

    ## Defining the self starter function  
    ssfct <- function(dframe)
    {   
        LL4values0 <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)
        LL4values <- LL4values0 * c(-1 / (mean(dframe[, 1])), 1, 1, 1)
        LL4values[notFixed]
    }

    ##Defining the first and second derivative (in the parameters)
    
    logGamma <- function(x, y)
    {
        lenX <- length(x)
        retVec <- rep(NA, lenX)
        for (i in 1:lenX)
        {
            intFct <- function(t){dgamma(t, y[i], 1) * log(t)}
            retVec[i] <- ifelse(x[i] < 1e-10, 0, integrate(intFct, 0, x[i])[[1]])
        }
        retVec
    }
            
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- pgamma(parmMat[, 1] * dose, parmMat[, 4], 1)

        cbind(
        t1 * dgamma(parmMat[, 1] * dose, parmMat[, 4], 1) * dose,
        1 - t2,
        t2,
        t1 * (logGamma(parmMat[, 1] * dose, parmMat[, 4]) - digamma(parmMat[, 4]) * t2)
        )[, notFixed]
    }
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        (parmMat[, 3] - parmMat[, 2]) * dgamma(parmMat[, 1] * x, parmMat[, 4], 1) * parmMat[, 1] 
    }

    ## Defining the ED function
    edfct <- NULL
  
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Gamma", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "gamma"
    invisible(returnList)
}
