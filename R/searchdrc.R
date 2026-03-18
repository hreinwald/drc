#' Search through a range of initial parameter values to obtain convergence
#'
#' \code{searchdrc} provides a facility for searching through a range of initial
#' values for a single parameter in order to obtain convergence of the non-linear
#' estimation procedure used in dose-response curve fitting.
#'
#' The function iterates through at most \code{len} evenly spaced values within
#' the specified \code{range}, using each as a starting value for the chosen
#' parameter. The search stops as soon as the first successful model fit is
#' found. You would need to identify the parameter which is most likely to cause
#' problems for the estimation procedure.
#'
#' Parameter names should be provided \strong{without} the curve suffix. For
#' example, use \code{"b"} rather than \code{"b:1"}. The function internally
#' matches the parameter using the pattern \code{"^<which>:"} against the full
#' parameter names stored in the model object.
#'
#' @param object an object of class \code{'drc'}, which must have valid
#'   \code{$start} and \code{$parNames} fields populated. This is typically
#'   an object from a model that failed to converge but was still constructed
#'   with initial parameter values.
#' @param which a character string containing the parameter name
#'   \strong{without} the curve suffix (e.g., \code{"b"} not \code{"b:1"}).
#'   Must exactly match one of the parameter names in the model object.
#' @param range a numeric vector of exactly length 2 specifying the interval
#'   endpoints \code{c(lower, upper)} for the search range. The two endpoints
#'   must be different.
#' @param len a positive integer (minimum 2). The maximum number of evenly
#'   spaced starting values to try within \code{range}. The search stops early
#'   as soon as convergence is achieved, so the actual number of attempts may
#'   be less than \code{len}. Defaults to \code{50}.
#' @param verbose logical. If \code{TRUE}, prints progress messages indicating
#'   which starting value is currently being tried. Defaults to \code{FALSE}.
#'
#' @return If convergence is achieved, returns the fitted model object of class
#'   \code{'drc'}, corresponding to the \strong{first} starting value in the
#'   search grid that led to a successful fit. If no starting value leads to
#'   convergence, the function throws an error.
#'
#' @author Christian Ritz, Hannes Reinwald.
#'
#' @seealso
#'   \code{\link[drc]{drm}} for the main model fitting function,
#'   \code{\link[drc]{drmc}} for control arguments,
#'   \code{\link[stats]{update}} for the update method used internally.
#'
#' @examples
#' \dontrun{
#' library(drc)
#'
#' # Fit an initial model (which may fail to converge)
#' myModel <- drm(response ~ dose, data = myData, fct = LL.4())
#'
#' # Search over a range of starting values for the slope parameter "b"
#' myModelFixed <- searchdrc(myModel, which = "b", range = c(-5, 5), len = 100)
#'
#' # With progress messages enabled
#' myModelFixed <- searchdrc(myModel, which = "b", range = c(-5, 5),
#'                           len = 100, verbose = TRUE)
#' }
#'
#' @keywords models nonlinear
searchdrc <- function(object, which, range, len = 50, verbose = FALSE)
{
  # 1. Input validation
  if (!inherits(object, "drc")) {
    stop("'object' must be of class 'drc'.", call. = FALSE)
  }
  if (is.null(object$start) || is.null(object$parNames)) {
    stop(paste0(
      "'object' must have valid '$start' and '$parNames' fields. ",
      "Ensure the model object was constructed with initial parameter values."
    ), call. = FALSE)
  }
  if (!is.character(which) || length(which) != 1 || nchar(trimws(which)) == 0) {
    stop("'which' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.numeric(range) || length(range) != 2) {
    stop("'range' must be a numeric vector of exactly length 2.", call. = FALSE)
  }
  if (range[1] == range[2]) {
    stop("The two endpoints of 'range' must be different.", call. = FALSE)
  }
  if (!is.numeric(len) || length(len) != 1 || len < 2) {
    stop("'len' must be a single numeric value of at least 2.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }
  
  len <- as.integer(len)
  
  # 2. Identify the target parameter index
  sv       <- object$start
  parNames <- object$parNames[[2]]
  
  matchPattern <- paste0("^", gsub("([.\\^$*+?\\[\\]\\{\\}()|])", "\\\\\\1", which), ":")
  matchIndices <- seq_along(parNames)[regexpr(matchPattern, parNames) > 0]
  
  if (length(matchIndices) == 0) {
    stop(paste0(
      "No parameter matching '", which, "' was found. ",
      "Available parameter names (without curve suffix) are: ",
      paste(unique(sub(":.*$", "", parNames)), collapse = ", "), "."
    ), call. = FALSE)
  }
  if (length(matchIndices) > 1) {
    warning(paste0(
      "Multiple parameters matched '", which, "': ",
      paste(parNames[matchIndices], collapse = ", "), ". ",
      "Using the first match: '", parNames[matchIndices[1]], "'."
    ), call. = FALSE)
  }
  
  # 3. Search loop
  on.exit(options(warn = getOption("warn")), add = TRUE)
  searchGrid <- seq(range[1], range[2], length.out = len)
  modelFit   <- NULL
  
  for (i in seq_along(searchGrid))
  {
    sv[matchIndices[1]] <- searchGrid[i]
    
    if (verbose) {
      message(sprintf("[searchdrc] Attempt %d / %d : %s = %.6g",
                      i, len, which, searchGrid[i]))
    }
    
    modelFit <- tryCatch(
      withCallingHandlers(
        update(object, start = sv, control = drmc(noMessage = TRUE)),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    
    if (!is.null(modelFit)) {
      if (verbose) {
        message(sprintf(
          "[searchdrc] Convergence achieved at %s = %.6g (attempt %d / %d).",
          which, searchGrid[i], i, len
        ))
      }
      break
    }
  }
  
  # 4. Return result or stop with informative error
  if (!is.null(modelFit)) {
    return(modelFit)
  }
  
  stop(paste0(
    "Convergence failed. No starting value for parameter '", which,
    "' in the range [", range[1], ", ", range[2], "] across ",
    len, " attempt(s) led to a successful fit. ",
    "Consider expanding the range or increasing 'len'."
  ), call. = FALSE)
}