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
#' the ordinary log-logistic model: \code{b1}, \code{b2} are the slopes and \code{e1}, \code{e2} the
#' locations of the two phases; \code{c} (here \code{c1}) is the lower asymptote of the whole curve
#' and the upper asymptote is \code{d1 + d2}. Starting values from the self starter may need to be
#' replaced by user-supplied \code{start} values. No effective doses (\code{\link{ED}}) are available.
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
#' @examples
#' ## Simulated two-phase data: fit with drm() and plot (starting values supplied)
#' dose <- rep(10^seq(-1, 4, length.out = 12), each = 2)
#' parm <- c(b1 = 2, c1 = 0.5, d1 = 4, e1 = 3, b2 = 1.5, d2 = 4, e2 = 400)
#' mu <- twophase()$fct(dose, matrix(parm, length(dose), 7, byrow = TRUE))
#' dat <- data.frame(dose = dose, response = mu + 0.04 * sin(seq_along(dose) * 1.7))
#' m1 <- drm(response ~ dose, data = dat, fct = twophase(), start = parm)
#' plot(m1, type = "all", log = "x", main = "Two-phase model")
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
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        fixed1.4 <- fixed[1:4]
        fixed5.7 <- fixed[5:7]
        LL.4()$fct(dose, parmMat[, 1:4, drop = FALSE]) + LL.3()$fct(dose, parmMat[, 5:7, drop = FALSE])
    }

    ## Defining self starter function
    ssfct <- function(dframe)
    {  
        first4 <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
    
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
