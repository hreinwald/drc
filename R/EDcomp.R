#' @title Comparison of relative potencies between dose-response curves
#'
#' @description
#' Relative potencies (also called selectivity indices) for arbitrary doses are compared between
#' fitted dose-response curves.
#'
#' @param object an object of class 'drc'.
#' @param percVec a numeric vector of dosage values.
#' @param percMat a matrix with 2 columns providing the pairs of indices of \code{percVec} to be
#'   compared. By default all pairs are compared.
#' @param compMatch an optional character vector of names of assays to be compared. If not specified
#'   all comparisons are supplied.
#' @param od logical. If TRUE adjustment for over-dispersion is used. This argument only makes a
#'   difference for binomial data.
#' @param vcov. function providing the variance-covariance matrix. \code{\link{vcov}} is the default,
#'   but \code{sandwich} is also an option (for obtaining robust standard errors).
#' @param reverse logical. If TRUE the order of comparison of two curves is reversed.
#' @param interval character string specifying the type of confidence intervals to be supplied.
#'   The default is \code{"none"}. Use \code{"delta"} for asymptotics-based confidence intervals,
#'   \code{"fieller"} for confidence intervals based on Fieller's theorem, or \code{"fls"} for
#'   confidence intervals back-transformed from logarithm scale.
#' @param level numeric. The level for the confidence intervals. Default is 0.95.
#' @param reference character string. Is the upper limit or the control level the reference?
#' @param type character string specifying whether absolute or relative response levels are supplied.
#' @param display logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).
#' @param pool logical. If TRUE curves are pooled. Otherwise they are not. This argument only works
#'   for models with independently fitted curves as specified in \code{\link{drm}}.
#' @param logBase numeric. The base of the logarithm in case logarithm transformed dose values are used.
#' @param multcomp logical to switch on output for use with the package \pkg{multcomp}. Default is FALSE.
#' @param ... additional arguments passed to the function doing the calculations.
#'
#' @return An invisible matrix containing the estimates and the corresponding estimated standard
#'   errors and possibly lower and upper confidence limits. Or, alternatively, a list with elements
#'   that may be plugged directly into \code{parm} in the package \pkg{multcomp} (when \code{multcomp}
#'   is TRUE).
#'
#' @details
#' Fieller's theorem is incorporated using the formulas provided by Kotz and Johnson (1983) and
#' Finney (1978).
#'
#' For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper'
#' argument or the 'interval' argument specifying limits for the bisection method.
#'
#' @seealso \code{\link{ED.drc}} for calculating effective doses.
#'
#' @examples
#' spinach.LL.4 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
#'
#' EDcomp(spinach.LL.4, c(50, 50))
#' EDcomp(spinach.LL.4, c(10, 50))
#' EDcomp(spinach.LL.4, c(10, 50), reverse = TRUE)
#'
#' @author Christian Ritz
#' @keywords models nonlinear
"EDcomp" <-
function(object, percVec, percMat = NULL, compMatch = NULL, od = FALSE, vcov. = vcov, reverse = FALSE, 
interval = c("none", "delta", "fieller", "fls"), level = ifelse(!(interval == "none"), 0.95, NULL), 
reference = c("control", "upper"), type = c("relative", "absolute"),
display = TRUE, pool = TRUE, logBase = NULL, multcomp = FALSE, ...)
{
     ## Matching the argument 'method'
     interval <- match.arg(interval)
     reference <- match.arg(reference)
     type <- match.arg(type)     

    if ( (is.null(logBase)) && (interval == "fls") )
    {
        stop("Argument 'logBase' not specified for interval = 'fls'")
    }

    ## Checking contents of percVec vector ... should be numbers between 0 and 100
    if ( (type == "relative") && any(percVec<=0 | percVec>=100) ) 
    {
        stop("Percentages outside the interval [0, 100] not allowed")
    }

    if (missing(compMatch)) {matchNames <- FALSE} else {matchNames <- TRUE}

    lenPV <- length(percVec)

    ## Retrieving relevant quantities
    indexMat <- object$"indexMat"
    parmMat <- object$"parmMat"


    curveNames <- colnames(object$"parmMat")
    if (any(suppressWarnings(is.na(as.numeric(curveNames)))))
    {
        curveOrder <- order(curveNames)
    } else { # if names are numbers then skip re-ordering
        curveOrder <- 1:length(curveNames)
    }
    
    strParm0 <- curveNames[curveOrder]
    indexMat <- indexMat[, curveOrder, drop = FALSE]
    lenEB <- ncol(indexMat) 
    sifct <- createsifct(object$"fct"$"edfct", logBase, identical(interval, "fls"), indexMat, length(coef(object)))    
    
    parmMat <- parmMat[, curveOrder, drop = FALSE]
    
    strParm <- strParm0
    varMat <- vcov.(object)

    ## Calculating SI values
    numComp <- (lenPV*(lenPV-1)/2)*(lenEB * (lenEB - 1) / 2)
    matchVec <- rep(TRUE, numComp)
    rNames <- rep("", numComp)
    oriMat <- matrix(0, numComp, 2)    
    degfree <- df.residual(object)  
    rowIndex <- 1

    pairsMat <- combinations(lenEB, 2)# canonical "2" as pairs are considered
    if (is.null(percMat))
    {
        percMat <- combinations(lenPV, 2)  # canonical "2" as pairs are considered
    }
    if (reverse)
    {
        pairsMat <- pairsMat[, 2:1, drop = FALSE]
        percMat <- percMat[, 2:1, drop = FALSE]
    }

    appFct1 <- function(percVal)
    {
        apply(pairsMat, 1, siInner, pVec = percVec[percVal], compMatch = compMatch, object = object, indexMat = indexMat, parmMat = parmMat, 
        varMat = varMat, level = level, reference = reference, type = type, sifct = sifct, interval = interval, degfree = degfree, logBase = logBase)
    }
    SImat0 <- matrix(apply(percMat, 1, appFct1), nrow = nrow(pairsMat) * nrow(percMat), byrow = TRUE)
    SImat <- SImat0[, 1:4, drop = FALSE]
    dSImat <- SImat0[, 5:ncol(SImat0), drop = FALSE]

    appFct2 <- function(percVal)
    {
        apply(pairsMat, 1, 
        function(indPair, percVal) 
        {
            paste(strParm[indPair[1]], "/", strParm[indPair[2]], ":", percVec[percVal[1]], "/", percVec[percVal[2]], sep = "")
        }, percVal = percVal)
    }    
    rownames(SImat) <- apply(percMat, 1, appFct2) 

    appFct3 <- function(percVal)
    {
        apply(pairsMat, 1, 
        function(indPair, percVal) 
        {
            (is.null(compMatch) || all(c(strParm[indPair[1]], strParm[indPair[2]]) %in% compMatch))
        })
    }
    SImat <- SImat[as.vector(apply(percMat, 1, appFct3)), , drop = FALSE]
    
    if (!identical(interval, "none"))
    {
        SImat <- SImat[, -4, drop = FALSE]
        cNames <- c("Estimate", "Lower", "Upper")
    
    } else {
        cNames <- c("Estimate", "Std. Error", "t-value", "p-value")
    }    
    colnames(SImat) <- cNames
    
    ciLabel <- switch(interval,
    "delta" = "Delta method",
    "tfls" = "To and from log scale",
    "fls" = "From log scale",
    "fieller" = "Fieller")
    
    resPrint(SImat, "Estimated ratios of effect doses", interval, ciLabel, display = display)
    
##    invisible(SImat) 

    if(multcomp)
    {         
        SImat1 <- SImat[, 1]
        namesVec <- names(SImat1)
        SImat1VC <- dSImat %*% varMat %*% t(dSImat)
        colnames(SImat1VC) <- namesVec
        rownames(SImat1VC) <- namesVec
      
        invisible(list(multcomp = parm(SImat1, SImat1VC)))
        
    } else {
        invisible(SImat)       
    } 
}


#' @title Fieller's confidence interval
#' @keywords internal
"fieller" <-
function(mu, df, vcMat, level = 0.95, finney = FALSE, resVar)
{
    tper <- qt(1-0.5*(1-level), df)^2 
    
    if (!finney)
    {
        ## Based on the entry on Fieller's theorem 
        ##  in Encyclopedia of Statistical Sciences Vol. 3 (1983), p. 86 
        ##  essentially same formula as in Finney (see below)

        mup <- prod(mu)
    
        fVec0 <- mup - tper*vcMat[1,2]
        y2 <- mu[2]^2    
        fVec <- (fVec0)^2 - (mu[1]^2 - tper*vcMat[1,1])*(y2 - tper*vcMat[2,2])
    
        denom <- y2 - tper*vcMat[2,2]
        lowerL <- (fVec0 - sqrt(fVec))/denom
        upperL <- (fVec0 + sqrt(fVec))/denom
        
    } else {
    
        ## Using the formula 
        ##  in Finney: Statistical Method in Biological Assay p. 81 (3rd edition, 1978)
        ##  OOPS: uses the estimated residual variance
        fac <- sqrt(tper)*sqrt(resVar)/mu[2]
        g <- tper*vcMat[2,2]/(mu[2]^2)
        if (g >= 1) {stop("Fieller's theorem not useful!")} 
        ratio <- mu[1]/mu[2]
    
        v11 <- vcMat[1,1]/(resVar)
        v12 <- vcMat[1,2]/(resVar)
        v22 <- vcMat[2,2]/(resVar)
        innerBr <- g*(v11 - (v12^2)/v22)
        inBr <- v11 - 2*ratio*v12 + (ratio^2)*v22 - innerBr

        firstTerm <- ratio - g*vcMat[1,2]/vcMat[2,2]
        secondTerm <- fac*sqrt(inBr)
        denom <- 1 - g
        lowerL <- (firstTerm - secondTerm)/denom
        upperL <- (firstTerm + secondTerm)/denom
    }
    return(c(lowerL, upperL))
}

#' @title Split index vectors into shared and unique components
#' @keywords internal
"splitInd"  <- function(ind1, ind2)
{
    matchVec1 <- ind1 %in% ind2
    matchVec2 <- ind2 %in% ind1
    lmv1 <- sum(matchVec1)
    if (lmv1 > 0.01)
    {
        inCommon <- matrix( c( (1:length(ind1))[matchVec1], (1:length(ind2))[matchVec2], ind1[matchVec1]), lmv1, 3)
    } else {
        inCommon <- NULL
    }
    
    only1 <- matrix( c( (1:length(ind1))[!matchVec1], ind1[!matchVec1] ), sum(!matchVec1), 2)

    only2 <- matrix( c( (1:length(ind2))[!matchVec2], ind2[!matchVec2] ), sum(!matchVec2), 2)

    return(list(only1, only2, inCommon))
}

#' @title Create selectivity index function
#' @keywords internal
createsifct <- function(edfct, logBase = NULL, fls = FALSE, indexMat, lenCoef)
{
    if (is.null(edfct)) 
    {
        stop("SI values cannot be calculated")
    } else {
        
        if (!fls)
        {
            if (is.null(logBase))  # this clause has been updated October 12 2010
            {
                "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
                {
                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- rep(0, lenCoef)
                    ED1d[indexMat[, jInd]] <- ED1[[2]]        
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- rep(0, lenCoef)
                    ED2d[indexMat[, kInd]] <- ED2[[2]]        

                    SIpair <- ED1v / ED2v
                    SIder <- (ED1d - SIpair * ED2d) / ED2v

                    return(list(val = SIpair, der = SIder,
                    der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
                }
            } else {
        
                "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
                {
                    ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                    ED1v <- ED1[[1]]
                    ED1d <- rep(0, lenCoef)
                    ED1d[indexMat[, jInd]] <- ED1[[2]]        
        
                    ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                    ED2v <- ED2[[1]]
                    ED2d <- rep(0, lenCoef)
                    ED2d[indexMat[, kInd]] <- ED2[[2]]        

                    SIpair <- logBase^(ED1v - ED2v)
                    SIder <- SIpair * log(logBase) * (ED1d - ED2d)

                    return(list(val = SIpair, der = SIder,
                    der1 = (log(logBase)*logBase^ED1v)*ED1d, der2 = (log(logBase)*logBase^ED2v)*ED2d, 
                    valnum = logBase^ED1v, valden = logBase^ED2v))
                }        
            }
        } else {
            
            "sifct" <- function(parm1, parm2, pair, jInd, kInd, reference, type, ...)
            {
                ED1 <- edfct(parm1, pair[1], reference = reference, type = type, ...)
                ED1v <- ED1[[1]]
                ED1d <- rep(0, lenCoef)
                ED1d[indexMat[, jInd]] <- ED1[[2]]        
        
                ED2 <- edfct(parm2, pair[2], reference = reference, type = type, ...)
                ED2v <- ED2[[1]]
                ED2d <- rep(0, lenCoef)
                ED2d[indexMat[, kInd]] <- ED2[[2]]        

                SIpair <- ED1v - ED2v
                SIder <- ED1d - ED2d

                return(list(val = SIpair, der = SIder,
                der1 = ED1d, der2 = ED2d, valnum = ED1v, valden = ED2v))
            }
        }        
        return(sifct)
    }
}

