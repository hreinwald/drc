#' @title Extract fitted values from model
#'
#' @description
#' Extracts fitted values from an object of class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param ... additional arguments.
#'
#' @return Fitted values extracted from \code{object}.
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' plot(fitted(ryegrass.m1), residuals(ryegrass.m1))  # a residual plot
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"fitted.drc" <-
function(object, ...)
{
#    if (missing(...))
#    {
#        return(object$"predres"[, 1])
#    } else {
#        predict(object, ...)
#    }
    predict(object, ...)
##    return(object$"predres"[, 1])
}
