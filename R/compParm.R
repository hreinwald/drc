#' @title Comparison of parameters
#'
#' @description
#' Compare parameters from different assays, either by means of ratios or differences.
#'
#' @param object an object of class 'drc'.
#' @param strVal a name of parameter to compare.
#' @param operator a character. If equal to \code{"/"} (default) parameter ratios are compared.
#'   If equal to \code{"-"} parameter differences are compared.
#' @param vcov. function providing the variance-covariance matrix. \code{\link{vcov}} is the default,
#'   but \code{sandwich} is also an option (for obtaining robust standard errors).
#' @param od logical. If TRUE adjustment for over-dispersion is used.
#' @param pool logical. If TRUE curves are pooled. Otherwise they are not. This argument only works
#'   for models with independently fitted curves as specified in \code{\link{drm}}.
#' @param display logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).
#'
#' @return A matrix with columns containing the estimates, estimated standard errors, values of
#'   t-statistics and p-values for the null hypothesis that the ratio equals 1 or that the difference
#'   equals 0 (depending on the \code{operator} argument).
#'
#' @seealso \code{\link{ED.drc}} for calculating effective doses and \code{\link{EDcomp}} for
#'   comparing effective doses.
#'
#' @examples
#' spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach,
#' fct = LL.4(names = c("b", "lower", "upper", "ed50")))
#'
#' ## Calculating ratios of parameter estimates for "ed50"
#' compParm(spinach.m1, "ed50")
#'
#' ## Calculating differences between parameter estimates for "ed50"
#' compParm(spinach.m1, "ed50", "-")
#'
#' @author Christian Ritz
#' @keywords models nonlinear
"compParm" <-
function(object, strVal, operator = "/", vcov. = vcov, od = FALSE, pool = TRUE, display = TRUE)
{
#    if (inherits(object, "mixdrc")) {sep <- ".{1}"} else {sep <- ":{1}"}
    sep <- ":{1}"
    presentVec <- grep(paste("^", strVal, sep, sep = ""), object$"parNames"[[1]])  # strParm)           

    lenPV <- length(presentVec)
    if (lenPV < 2) 
    {
        stop("No parameters to compare")
    }

    ## Extracting information from model fit 
#    if (inherits(object, "mixdrc")) 
#    {
#        sumObj <- summary(object)
#        parm <- sumObj$"coefficients"
#        varMat <- sumObj$"varMat"
#    } else {
#        parm <- as.vector(coef(object))
#        varMat <- vcov(object, od = od, pool = pool)
#    }
    parm <- as.vector(coef(object))
    varMat <- vcov.(object)
    
    ## Defining comparison function and its derivative
    if (identical(operator, "/"))
    {
        hypVal <- 1
        fct <- function(ind){parm[ind[1]] / parm[ind[2]]}
        dfct <- function(ind)
        {
            transVec <- c(1 / parm[ind[2]], -parm[ind[1]] / (parm[ind[2]]^2))
            sqrt(transVec %*% varMat[ind,ind] %*% transVec)
        }
    }
    if (identical(operator, "-"))
    {
        hypVal <- 0
        fct <- function(ind){parm[ind[1]] - parm[ind[2]]}
        transVec <- c(1, -1)
        dfct <- function(ind){sqrt(transVec %*% varMat[ind,ind] %*% transVec)}
    }

    ## Calculating differences or ratios
    lenRV <- lenPV*(lenPV-1)/2
    cpMat <- matrix(0, lenRV, 4)
    compParm <- rep("", lenRV)
    
    degfree <- df.residual(object)
    if (is.null(degfree)) {degfree <- 100}  # ad hoc solution for mixdrc    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (object$"type" == "continuous")
    {
        pFct <- function(x) {pt(x, degfree)}
    } else {
        pFct <- pnorm
    }        
   
    strParm <- object$"parNames"[[3]]
    k <- 1
    for (i in 1:lenPV) 
    {
        for (j in 1:lenPV)
        {
            if (j<=i) {next}

            cpMat[k, 1] <- fct(presentVec[c(i,j)])  #parm[i]/parm[j]
            cpMat[k, 2] <- dfct(presentVec[c(i,j)])

            tVal <- (cpMat[k, 1] - hypVal)/cpMat[k, 2]
            cpMat[k, 3] <- tVal
            cpMat[k, 4] <- pFct(-abs(tVal)) + (1 - pFct(abs(tVal)))

            compParm[k] <- paste(strParm[presentVec[c(i, j)]], collapse = operator)
            k <- k+1
        }
    }
    dimnames(cpMat) <- list(compParm, c("Estimate", "Std. Error", "t-value", "p-value"))
    
    if (display)
    {
        cat("\nComparison of parameter", paste("'", strVal, "'", sep = ""), "\n\n")
        printCoefmat(cpMat)
    }
    invisible(cpMat)
}
