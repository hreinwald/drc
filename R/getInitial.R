#' Showing starting values used
#'
#' Returns the starting values of the model parameters used when fitting a dose-response model.
#'
#' @param object object of class 'drc'.
#'
#' @return A vector of starting values for the model parameters used to initialize the
#'   estimation procedure.
#'
#' @author Christian Ritz
#'
#' @note This function is masking the standard function in the stats package.
#'
#' @keywords models nonlinear
"getInitial" <- function(object)
{
    initval <- object$"start"
    names(initval) <- object$"parNames"[[2]]
    
    initval
}