#' Convert absolute to relative response levels
#'
#' Internal helper that converts an absolute response level to a relative (percentage) scale
#' based on the upper and lower asymptotes of a dose-response curve.
#'
#' @param parmVec numeric vector of model parameters where the third element is the upper
#'   asymptote and the second element is the lower asymptote.
#' @param respl numeric response level to convert.
#' @param typeCalc character string. If "absolute", the conversion is performed;
#'   otherwise the input \code{respl} is returned unchanged.
#'
#' @return A numeric value representing the (possibly converted) response level as a percentage.
#'
#' @keywords internal
"absToRel" <- function(parmVec, respl, typeCalc)
{
    ## Converting absolute to relative
    if (typeCalc == "absolute") 
    {
        denom <- parmVec[3] - parmVec[2]
        if (abs(denom) < .Machine$double.eps)
        {
            stop("Cannot convert absolute to relative response: upper and lower asymptotes are equal")
        }
        p <- 100 * ((parmVec[3] - respl) / denom)
    } else {  
        p <- respl
    }

    p
}