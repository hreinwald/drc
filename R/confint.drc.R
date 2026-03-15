#' @title Confidence Intervals for Model Parameters
#'
#' @description
#' Computes confidence intervals for one or more parameters in a fitted
#' dose-response model of class `"drc"`. Confidence intervals are constructed
#' using either a t-distribution (for continuous response models) or a standard
#' normal distribution (for all other response types).
#'
#' @param object A fitted model object of class `"drc"`.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of indices or a vector of parameter name strings.
#'   If missing, all parameters are considered.
#' @param level The confidence level required. Defaults to `0.95`.
#' @param pool Logical. If `TRUE` (default), curves are pooled. Otherwise they
#'   are not. This argument only works for models with independently fitted
#'   curves as specified in [drm()].
#' @param ... Additional arguments for methods. Currently not used.
#'
#' @return A numeric matrix with two columns giving the lower and upper
#'   confidence limits for each parameter. Columns are labelled as
#'   \eqn{\frac{(1 - \text{level})}{2} \times 100\%} and
#'   \eqn{\left(1 - \frac{(1 - \text{level})}{2}\right) \times 100\%}
#'   (by default \code{2.5 \%} and \code{97.5 \%}).
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [drm()] — for fitting dose-response models.
#' * [confint.basic()] — the internal helper used to construct the intervals.
#' * [summary.drc()] — for a full summary of model coefficients.
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' ## Confidence intervals for all parameters
#' confint(ryegrass.m1)
#'
#' ## Confidence interval for a single parameter
#' confint(ryegrass.m1, "e")
#'
#' @keywords models nonlinear
#' @export
"confint.drc" <- function(object, parm, level = 0.95, pool = TRUE, ...)
{
  ## Matching parameter names
  if (!missing(parm))
  {
    matchVec <- object$"parNames"[[2]] %in% parm
    if (!any(matchVec))
    {
      stop("The 'parm' argument does not match an actual parameter name.")
    }
  } else {
    matchVec <- rep(TRUE, length(object$"parNames"[[2]]))
  }
  
  ## Constructing matrix of confidence intervals
  confint.basic(
    summary(object, pool = pool)$"coefficients"[matchVec, 1:2, drop = FALSE],
    level,
    object$"type",
    df.residual(object)
  )
}


#' @title Basic Confidence Interval Calculation
#'
#' @description
#' An internal helper function that constructs a confidence interval matrix
#' from a matrix of parameter estimates and their standard errors. A
#' t-distribution quantile is used for continuous response models; a standard
#' normal quantile is used for all other response types (binomial, event,
#' Poisson, negbin1, negbin2).
#'
#' @param estMat A numeric matrix with two columns: the first column contains
#'   parameter estimates and the second column contains their standard errors.
#' @param level The confidence level required (e.g., `0.95` for 95% intervals).
#' @param intType A character string specifying the response type of the model.
#'   One of `"binomial"`, `"continuous"`, `"event"`, `"Poisson"`,
#'   `"negbin1"`, or `"negbin2"`. Determines whether a normal or t-distribution
#'   quantile is used. For `"continuous"` models a t-distribution with `dfres`
#'   degrees of freedom is used; all other types use the standard normal.
#' @param dfres The residual degrees of freedom. Only used when
#'   `intType = "continuous"`.
#' @param formatting Logical. If `TRUE` (default), row and column names are
#'   added to the returned matrix.
#'
#' @return A numeric matrix with two columns giving the lower and upper
#'   confidence limits for each parameter.
#'
#' @seealso [confint.drc()] — the user-facing function that calls this helper.
#'
#' @keywords internal
"confint.basic" <- function(estMat, level, intType, dfres, formatting = TRUE)
{
  alphah <- (1 - level) / 2
  
  tailPercentile <- switch(intType,
                           binomial   = qnorm(1 - alphah),
                           continuous = qt(1 - alphah, dfres),
                           event      = qnorm(1 - alphah),
                           Poisson    = qnorm(1 - alphah),
                           negbin1    = qnorm(1 - alphah),
                           negbin2    = qnorm(1 - alphah),
                           stop(paste0(
                             "Unknown intType '", intType, "'. ",
                             "Must be one of: 'binomial', 'continuous', 'event', ",
                             "'Poisson', 'negbin1', 'negbin2'."
                           ))
  )
  
  estVec     <- estMat[, 1]
  halfLength <- tailPercentile * estMat[, 2]
  confMat    <- matrix(c(estVec - halfLength, estVec + halfLength), ncol = 2)
  
  if (formatting)
  {
    colnames(confMat) <- c(
      paste(format(100 * alphah),       "%", sep = " "),
      paste(format(100 * (1 - alphah)), "%", sep = " ")
    )
    rownames(confMat) <- rownames(estMat)
  }
  
  return(confMat)
}