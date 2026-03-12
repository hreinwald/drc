#' Maximum mean response
#'
#' Estimates the maximum mean response and the dose at which it occurs. This function is only
#' implemented for the built-in functions of class \code{\link{braincousens}} and \code{\link{cedergreen}}.
#'
#' @param object an object of class 'drc'.
#' @param lower numeric. Lower limit for bisection method. Must be smaller than the EDx level to be calculated.
#' @param upper numeric. Upper limit for bisection method. Must be larger than the EDx level to be calculated.
#' @param pool logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for
#'   models with independently fitted curves as specified in \code{\link{drm}}.
#'
#' @return A matrix with one row per curve in the data set and two columns: one containing the dose
#'   at which the maximum occurs and one containing the corresponding maximum response.
#'
#' @references Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved empirical models
#'   describing hormesis, \emph{Environmental Toxicology and Chemistry} \bold{24}, 3166--3172.
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Fitting a Cedergreen-Ritz-Streibig model
#' lettuce.m1 <- drm(weight~conc, data = lettuce, fct = CRS.4c())
#'
#' ## Finding maximum average response and the corresponding dose
#' MAX(lettuce.m1)
#'
#' @keywords models nonlinear
"MAX" <- function(
object, lower = 1e-3, upper = 1000, pool = TRUE)
{
    ## Checking class of 'object'
    MAXlist <- object[[11]]$"maxfct"
    if (is.null(MAXlist)) {stop("No method available")}

    ## Retrieving relevant quantities
    indexMat <- object$"indexMat"
    parm <- as.vector(coef(object))
    parmMat <- object$"parmMat"
    strParm <- colnames(parmMat)    
    varMat <- vcov(object, pool = pool)

    ## Calculating ED values    
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM    
    dimNames <- rep("", ncolIM)
    MAXmat <- matrix(0, ncolIM, 2)

    for (i in indexVec)
    {
        parmInd <- indexMat[,i]
        varCov <- varMat[parmInd, parmInd]
        parmChosen <- parmMat[,i]
        MAXmat[i, ] <- MAXlist(parmChosen, lower, upper)
        dimNames[i] <- strParm[i]
    }

    dimnames(MAXmat) <- list(dimNames, c("Dose", "Response"))
    printCoefmat(MAXmat)
    invisible(MAXmat)    
}
