#' Calculation of backfit values from a fitted dose-response model
#'
#' By inverse regression backfitted dose values are calculated for the mean response per dose.
#'
#' @param drcObject an object of class 'drc'.
#'
#' @return Two columns with the original dose values and the corresponding backfitted values
#'   using the fitted dose-response model. For extreme dose values (e.g., high dose) the
#'   backfitted values may not be well-defined.
#'
#' @author Christian Ritz after a suggestion from Keld Sorensen.
#'
#' @seealso A related function is \code{\link{ED.drc}}.
#'
#' @examples
#' ryegrass.LL.4 <- drm(rootl~conc, data=ryegrass, fct=LL.4())
#'
#' backfit(ryegrass.LL.4)
#'
#' @keywords models nonlinear
backfit <- function(drcObject)
{
    DL <- drcObject$dataList
    DLdose <- DL$dose
    meansVec <- tapply(DL$origResp, DLdose, mean, na.rm = TRUE) 
    # arranged according to ascending dose values
    # therefore unique doses are sorted below

    backfitValues <- ED(drcObject, meansVec, type = "absolute", 
                        display = FALSE, multcomp = FALSE)[, 1, drop = FALSE]

#     colnames(backfitValues) <- "backfit"
#     rownames(backfitValues) <- sort(unique(DLdose))
#     backfitValues
    
    retMat <- cbind(dose = sort(unique(DLdose)), backfit = backfitValues)
    rownames(retMat) <- NULL
    return(retMat)
}