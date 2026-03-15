#' R-squared for dose-response models
#'
#' Calculates and displays R-squared values for a fitted dose-response model. For models
#' with multiple curves, per-curve and total R-squared values are returned.
#'
#' @param object an object of class 'drc'.
#'
#' @return Invisibly returns a matrix of R-squared values. For single-curve models, a 1x1 matrix.
#'   For multi-curve models, includes per-curve values and a total R-squared.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"Rsq" <- function(object)
{
    response <- object$data[,2]
    curve <- object$data[,4]
    uniCurve <- unique(curve)
    lenUC <- length(uniCurve)

    numerator <- tapply( residuals(object)^2, curve, sum)  # residual ss
    denominator <- tapply( (response - mean(response))^2, curve, sum)  # total SS

    totnum <- sum(residuals(object)^2)
    totden <- sum((response - mean(response))^2)

    ## Handle zero denominator (constant response) to avoid NaN
    rsqVals <- ifelse(denominator == 0, NA_real_, 1 - numerator / denominator)
    totRsq <- ifelse(totden == 0, NA_real_, 1 - totnum / totden)

    if (lenUC==1)
    {
        hText <- "\nR-square value\n"
        rsq <- matrix(rsqVals, 1, 1)
        rownames(rsq) <- "" 
    } else {
        hText <- "\nR-square values\n"
        rsq <- matrix(c(rsqVals, totRsq), lenUC+1, 1)
        rownames(rsq) <- c(as.character(uniCurve), "Total") 
    }
    colnames(rsq) <- ""
    
    cat(hText)
    printCoefmat(rsq)
    invisible(rsq)
}
