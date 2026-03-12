#' Expected or predicted response
#'
#' Returns the expected or predicted response for specified dose values. This is a
#' convenience function for easy access to predicted values.
#'
#' @param object object of class \code{drc} obtained from fitting a dose-response model.
#' @param xVec numeric vector of dose values.
#' @param ... additional arguments passed to \code{\link[drc]{predict.drc}}.
#'
#' @return A numeric vector of predicted values or possibly a matrix of predicted values
#'   and corresponding standard errors.
#'
#' @author Christian Ritz after a suggestion from Andrew Kniss.
#'
#' @seealso \code{\link[drc]{predict.drc}}
#'
#' @examples
#' ryegrass.m1 <- drm(ryegrass, fct = LL.4())
#' PR(ryegrass.m1, c(5, 10))
#'
#' @keywords models nonlinear
"PR" <- function(object, xVec, ...)
{
    lenXV <- length(xVec)
    
    curveId <- as.character(unique(object$data[, 3]))
    lenCI <- length(curveId)
    
    if (lenCI > 1)
    {
        retMat <- predict(object, data.frame(xVec, rep(curveId, rep(lenXV, lenCI))), se.fit = TRUE, ...)
        rownames(retMat) <- paste(rep(curveId, rep(lenXV, lenCI)), rep(as.character(xVec), lenCI), sep = ":")
    } else {
        retMat <- predict(object, data.frame(xVec))
        if (is.matrix(retMat))
        {
            rownames(retMat) <- rep(as.character(xVec), lenCI)
        } else {
            names(retMat) <- rep(as.character(xVec), lenCI)
        }
    }
    
    return(retMat)
}
