#' Residual sum of squares for dose-response models
#'
#' Calculates and displays the residual sum of squares (RSS) for a fitted dose-response model.
#' For models with multiple curves, per-curve and total RSS values are returned.
#'
#' @param object an object of class 'drc'.
#'
#' @return Invisibly returns a matrix of RSS values. For single-curve models, a 1x1 matrix.
#'   For multi-curve models, includes per-curve values and a total RSS.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
#' @export
"rss" <- function(object)
{
    curve <- object$data[,4]
    uniCurve <- unique(curve)
    lenUC <- length(uniCurve)

    rssVals <- tapply(residuals(object)^2, curve, sum)
    totRSS <- sum(residuals(object)^2)

    if (lenUC == 1)
    {
        hText <- "\nResidual sum of squares\n"
        rssMat <- matrix(rssVals, 1, 1)
        rownames(rssMat) <- ""
    } else {
        hText <- "\nResidual sums of squares\n"
        rssMat <- matrix(c(rssVals, totRSS), lenUC + 1, 1)
        rownames(rssMat) <- c(as.character(uniCurve), "Total")
    }
    colnames(rssMat) <- ""

    cat(hText)
    printCoefmat(rssMat)
    invisible(rssMat)
}
