#' @title Extract Model Coefficients
#'
#' @description
#' Extract parameter estimates.
#'
#' @param object an object of class 'drc'.
#' @param ... additional arguments.
#'
#' @return A vector of parameter coefficients which are extracted from the
#'   model object \code{object}.
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' coef(ryegrass.m1)
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"coef.drc" <-
function(object, ...)
{
    if (!is.null(object$"coefficients"))
    {
        return(object$"coefficients")
    } else {
        retVec <- object$fit$par
        names(retVec) <- object$parNames[[1]]
        return(retVec)
    }
}
