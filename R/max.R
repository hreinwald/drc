#' Maximum mean response
#'
#' Estimates the maximum mean response and the dose at which it occurs, using a
#' bisection method to locate the peak of the fitted dose-response curve. This
#' function is only implemented for the built-in model functions of class
#' \code{\link{braincousens}} and \code{\link{cedergreen}}, which are capable of
#' exhibiting hormesis (i.e., a non-monotone response with a stimulatory effect
#' at low doses).
#'
#' @param object an object of class \code{drc}, fitted using \code{\link{drm}}
#'   with a hormesis model such as \code{\link{CRS.4c}} or \code{\link{BC.4}}.
#'
#' @param lower numeric. Lower bound of the interval used by the bisection
#'   method to search for the dose at maximum response. Must be strictly smaller
#'   than \code{upper} and should be set below the expected dose at maximum
#'   response. Defaults to \code{1e-3}.
#'
#' @param upper numeric. Upper bound of the interval used by the bisection
#'   method to search for the dose at maximum response. Must be strictly larger
#'   than \code{lower} and should be set above the expected dose at maximum
#'   response. Defaults to \code{1000}.
#'
#' @param pool logical. If \code{TRUE} (default), curves are pooled when
#'   computing the variance-covariance matrix. Otherwise they are not. This
#'   argument only works for models with independently fitted curves as
#'   specified in \code{\link{drm}}. Note: currently the variance-covariance
#'   matrix is retrieved for internal consistency but standard errors are not
#'   yet reported in the output.
#'
#' @return Invisibly returns a numeric matrix with one row per curve in the
#'   data set and two columns:
#'   \describe{
#'     \item{Dose}{The dose at which the maximum mean response occurs, found
#'       via bisection within \code{[lower, upper]}.}
#'     \item{Response}{The estimated maximum mean response at that dose.}
#'   }
#'   Row names correspond to curve identifiers. If the computation fails for a
#'   given curve, the corresponding row will contain \code{NA} values and a
#'   warning is issued. The matrix is also printed to the console via
#'   \code{\link{printCoefmat}}.
#'
#' @details
#' The function numerically locates the dose \eqn{d^*} that maximises the fitted
#' dose-response curve over the search interval \code{[lower, upper]}:
#' \deqn{d^* = \arg\max_{d} f(d, \hat{\theta})}
#' where \eqn{f} is the fitted dose-response function and \eqn{\hat{\theta}} is
#' the vector of estimated parameters. The search is performed using a bisection
#' approach defined internally by the model's \code{maxfct} component.
#'
#' It is the user's responsibility to ensure that the true maximum lies within
#' \code{[lower, upper]}. If the maximum falls outside this interval, the
#' function will silently return a boundary value and a warning is issued.
#'
#' @references
#' Cedergreen, N., Ritz, C., and Streibig, J. C. (2005) Improved empirical
#' models describing hormesis, \emph{Environmental Toxicology and Chemistry}
#' \bold{24}, 3166--3172.
#'
#' @author Christian Ritz. Issues fixed and documentation enhanced by Hannes Reinwald.
#'
#' @examples
#' ## Fitting a Cedergreen-Ritz-Streibig model
#' lettuce.m1 <- drm(weight ~ conc, data = lettuce, fct = CRS.4c())
#'
#' ## Finding the maximum mean response and the corresponding dose
#' MAX(lettuce.m1)
#'
#' ## Custom search interval
#' MAX(lettuce.m1, lower = 1e-5, upper = 500)
#'
#' ## Capture the result matrix
#' result <- MAX(lettuce.m1)
#' result["Dose"]
#'
#' @keywords models nonlinear
"MAX" <- function(object, lower = 1e-3, upper = 1000, pool = TRUE)
{

  #  1. Validate class of 'object'
  if (!inherits(object, "drc")) {
    stop("'object' must be of class 'drc'. Please supply a model fitted with drm().")
  }

  
  #  2. Check that the model supports MAX (braincousens / cedergreen)
  MAXlist <- object[["fct"]][["maxfct"]]
  if (is.null(MAXlist)) {
    stop(paste(
      "No 'maxfct' method available for this model.",
      "MAX() is only supported for 'braincousens' and 'cedergreen' model classes."
    ))
  }

  
  #  3. Validate lower / upper bounds
  if (!is.numeric(lower) || length(lower) != 1 || !is.finite(lower)) {
    stop("'lower' must be a single finite numeric value.")
  }
  if (!is.numeric(upper) || length(upper) != 1 || !is.finite(upper)) {
    stop("'upper' must be a single finite numeric value.")
  }
  if (lower >= upper) {
    stop(paste0(
      "'lower' (", lower, ") must be strictly less than 'upper' (", upper, ")."
    ))
  }

  
  #  4. Retrieve relevant model components
  indexMat <- object[["indexMat"]]
  parmMat  <- object[["parmMat"]]
  strParm  <- colnames(parmMat)
  
  # vcov is retrieved for future SE support; pool affects its computation
  varMat <- tryCatch(
    vcov(object, pool = pool),
    error = function(e) {
      warning("Could not compute variance-covariance matrix: ", e$message)
      NULL
    }
  )
  

  #  5. Initialise output — use NA to distinguish failure from zero
  ncolIM   <- ncol(indexMat)
  indexVec <- seq_len(ncolIM)
  dimNames <- vector("character", ncolIM)
  MAXmat   <- matrix(NA_real_, nrow = ncolIM, ncol = 2)
  

  #  6. Loop over each curve
  for (i in indexVec)
  {
    # Curve label — fall back gracefully if missing
    curveName    <- if (!is.null(strParm) && !is.na(strParm[i])) strParm[i] else paste0("Curve_", i)
    dimNames[i]  <- curveName
    parmChosen   <- parmMat[, i]
    
    # Attempt to compute maximum for this curve
    MAXmat[i, ] <- tryCatch(
      {
        result <- MAXlist(parmChosen, lower, upper)
        
        # Warn if the optimum landed on a boundary (likely out-of-range)
        # Use unname() so that named return values (e.g. from cedergreen) are
        # compared correctly with the unnamed lower/upper scalars.
        # A tolerance of 1e-3 is used because numerical optimisers (e.g.
        # optimize()) return values near—but not exactly at—the boundary.
        bnd_tol <- 1e-3
        if (isTRUE(all.equal(unname(result[1]), lower, tolerance = bnd_tol)) ||
            isTRUE(all.equal(unname(result[1]), upper, tolerance = bnd_tol))) {
          warning(
            "The estimated maximum dose for curve '", curveName,
            "' is at the boundary of [lower, upper] = [", lower, ", ", upper, "]. ",
            "Consider widening the search interval."
          )
        }
        result
      },
      error = function(e) {
        warning(
          "MAX computation failed for curve '", curveName, "': ", e$message,
          ". Returning NA for this curve."
        )
        c(NA_real_, NA_real_)
      }
    )
  }

  
  #  7. Format and return
  dimnames(MAXmat) <- list(dimNames, c("Dose", "Response"))
  printCoefmat(MAXmat, na.print = "NA")
  invisible(MAXmat)
}
