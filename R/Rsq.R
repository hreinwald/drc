#' R-squared for dose-response models
#'
#' Calculates and displays R-squared values for a fitted dose-response model. For models
#' with multiple curves, per-curve and total R-squared values are returned.
#'
#' R-squared is computed as \eqn{1 - RSS / TSS} where RSS is the residual sum of squares
#' (obtained via [rss()]) and TSS is the total sum of squares.
#'
#' @param object an object of class 'drc'.
#'
#' @return Invisibly returns a matrix of R-squared values. For single-curve models, a 1x1 matrix.
#'   For multi-curve models, includes per-curve values and a total R-squared.
#'
#' @seealso [rss()] for the underlying residual sum of squares.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
#' @export
"Rsq" <- function(object)
{
    response <- object$data[,2]
    curve <- object$data[,4]
    uniCurve <- unique(curve)
    lenUC <- length(uniCurve)

    ## Use rss() for the residual sum of squares
    rssMat <- rss(object, print = FALSE)
    numerator <- rssMat[seq_len(lenUC), 1]
    denominator <- tapply( (response - mean(response))^2, curve, sum)  # total SS

    totnum <- if (lenUC == 1) numerator else rssMat[lenUC + 1, 1]
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
