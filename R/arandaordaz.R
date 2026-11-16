#' Asymptotic Regression Model
#'
#' The base function for the asymptotic regression model, providing the mean
#' function and self starter for a three-parameter model.
#'
#' The asymptotic regression model is a three-parameter model with mean function:
#'
#' \deqn{f(x) = a + (b-a)(1-\exp(-c x))}
#'
#' The parameter \eqn{a} is the lower limit (the response at \eqn{x=0}), \eqn{b} is the upper limit
#' (asymptote) and \eqn{c>0} is the rate of the increase. The default \code{names} are
#' \code{c("a", "b", "c")}. The equivalent models \code{\link{AR.3}} and \code{\link{AR.2}}
#' (parameterised by the scale \eqn{e = 1/c}) are usually preferred.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what
#'   value they are fixed. Use \code{NA} for parameters that are not fixed.
#'   Must be of length 3 (order: lower limit \code{a}, upper limit \code{b}, rate \code{c}).
#' @param names character vector of length 3 giving the names of the parameters
#'   (should not contain ":").
#' @param fctName optional character string used internally by convenience
#'   functions. Defaults to \code{"arandaordaz"} if not provided.
#' @param fctText optional character string used internally by convenience
#'   functions. Defaults to \code{"Asymptotic regression"} if not provided.
#'
#' @return A list of class \code{drcMean} with the following components:
#'   \describe{
#'     \item{fct}{The mean function taking arguments \code{dose} and \code{parm}.}
#'     \item{ssfct}{Self-starter function for generating initial parameter
#'       estimates from data.}
#'     \item{names}{Character vector of non-fixed parameter names.}
#'     \item{deriv1}{Reserved first derivative slot (currently \code{NULL}).}
#'     \item{deriv2}{Reserved second derivative slot (currently \code{NULL}).}
#'     \item{derivx}{Reserved derivative-with-respect-to-x slot (currently
#'       \code{NULL}).}
#'     \item{edfct}{Function for calculating effective dose (ED) values and
#'       their derivatives.}
#'     \item{inversion}{Inverse mean function for back-calculating dose from
#'       response.}
#'     \item{name}{Character string identifying the model function name.}
#'     \item{text}{Character string with a human-readable model description.}
#'     \item{noParm}{Integer giving the number of non-fixed parameters.}
#'   }
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso \code{\link{AR.2}}, \code{\link{AR.3}}, \code{\link{EXD.2}},
#'   \code{\link{EXD.3}}
#'
#' @examples
#' ## Simulated increasing data: fit with drm() and plot
#' dose <- rep(c(0.5, 1, 2, 4, 8, 16, 32), each = 3)
#' mu <- arandaordaz()$fct(dose, matrix(c(1, 8, 0.15), length(dose), 3, byrow = TRUE))
#' dat <- data.frame(dose = dose, response = mu + 0.15 * sin(seq_along(dose) * 2.3))
#' m1 <- drm(response ~ dose, data = dat, fct = arandaordaz())
#' plot(m1, type = "all", main = "Asymptotic regression")
#'
#' ## Effect of the lower limit a (fixed at 0) and the upper asymptote b (fixed at 6)
#' m.a0 <- drm(response ~ dose, data = dat, fct = arandaordaz(fixed = c(0, NA, NA)))
#' m.b6 <- drm(response ~ dose, data = dat, fct = arandaordaz(fixed = c(NA, 6, NA)))
#' plot(m1, type = "all", main = "Fixed lower limit and asymptote")
#' plot(m.a0, type = "none", add = TRUE, col = "darkgreen")
#' plot(m.b6, type = "none", add = TRUE, col = "red")
#' legend("bottomright", c("free", "a = 0", "b = 6"), col = c("black", "darkgreen", "red"),
#'   lty = 1, bty = "n")
#'
#' @keywords models nonlinear
#' @export
arandaordaz <- function(
    fixed   = c(NA, NA, NA),
    names   = c("a", "b", "c"),
    fctName,
    fctText
) {
  
  ## --- Input validation
  numParm <- 3

  if (length(fixed) != numParm) {
    stop("'fixed' must have length ", numParm)
  }
  if (is.list(fixed) || (!is.numeric(fixed) && !all(is.na(fixed)))) {
    stop("'fixed' must be a numeric vector")
  }
  if (!is.character(names) || length(names) != numParm) {
    stop("'names' must be a character vector of length ", numParm)
  }

  ## --- Handling 'fixed' argument
  # Convert to numeric if all NA (default c(NA, NA, NA) is logical)
  if (all(is.na(fixed))) {
    fixed <- as.numeric(fixed)
  }
  notFixed <- is.na(fixed)
  parmVec  <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]
  
  ## --- Mean function
  fct <- function(dose, parm) {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    
    parmMat[, 1] + (parmMat[, 2] - parmMat[, 1]) * (1 - exp(-parmMat[, 3] * dose))
  }
  
  ## --- Self-starter function
  # Nudge factors keep initial bounds just outside the observed data range,
  # preventing boundary issues in log transformations during starting-value
  # estimation.
  LOWER_SHRINK <- 0.95
  UPPER_EXPAND <- 1.05
  
  ssfct <- function(dataf) {
    x <- dataf[, 1]
    y <- dataf[, 2]
    
    aPar <- min(y) * LOWER_SHRINK
    bPar <- max(y) * UPPER_EXPAND
    
    # Compute pseudo-y values for log-linearisation; guard against
    # non-positive values that would produce NaN or -Inf from log().
    innerVal <- -((y - aPar) / (bPar - aPar) - 1)
    if (any(innerVal <= 0)) {
      warning(
        "Self-starter encountered invalid log argument; ",
        "initial estimates may be unreliable."
      )
      innerVal <- pmax(innerVal, .Machine$double.eps)
    }
    
    pseudoY <- log(innerVal)
    
    # Linear regression through the origin on pseudo-y values to
    # estimate the rate parameter.
    cPar <- coef(lm(pseudoY ~ I(-x) - 1))
    
    c(aPar, bPar, cPar)[notFixed]
  }
  
  ## --- Subset names to non-fixed parameters
  names <- names[notFixed]
  
  ## --- Effective dose function
  edfct <- function(parm, respl, reference, type, ...) {
    
    # Build the full parameter vector from estimated and fixed components.
    localParmVec <- parmVec
    localParmVec[notFixed] <- parm
    
    # Convert to relative scale if the type is absolute.
    if (type == "absolute") {
      p <- 100 * ((localParmVec[2] - respl) / (localParmVec[2] - localParmVec[1]))
    } else {
      p <- respl
    }
    
    # Adjust reference direction if calculated relative to the control.
    if (reference == "control") {
      p <- 100 - p
    }
    
    pProp <- p / 100
    EDp   <- -log(pProp) / localParmVec[3]
    
    EDder <- c(0, 0, log(pProp) / (localParmVec[3]^2))
    
    list(EDp, EDder[notFixed])
  }
  
  ## --- Inverse function 
  invfct <- function(y, parm) {
    localParmVec <- parmVec
    localParmVec[notFixed] <- parm
    
    log(-((y - localParmVec[1]) / (localParmVec[2] - localParmVec[1]) - 1)) /
      (-localParmVec[3])
  }
  
  ## --- Assemble and return the model object 
  returnList <- list(
    fct       = fct,
    ssfct     = ssfct,
    names     = names,
    deriv1    = NULL,
    deriv2    = NULL,
    derivx    = NULL,
    edfct     = edfct,
    inversion = invfct,
    name      = if (missing(fctName)) "arandaordaz" else fctName,
    text      = if (missing(fctText)) "Asymptotic regression" else fctText,
    noParm    = sum(is.na(fixed))
  )
  
  class(returnList) <- "drcMean"
  invisible(returnList)
}
