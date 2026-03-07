#' Estimation of ED values using model-averaging
#'
#' Estimates and confidence intervals for ED values are estimated using model-averaging.
#'
#' Model-averaging of individual estimates is carried out as described by Buckland
#' \emph{et al.} (1997) and Kang \emph{et al.} (2000) using AIC-based weights. The two
#' approaches differ w.r.t. the calculation of confidence intervals: Buckland \emph{et al.}
#' (1997) provide an approximate variance formula under the assumption of perfectly correlated
#' estimates (so, confidence intervals will tend to be too wide). Kang \emph{et al.} (2000)
#' use the model weights to calculate confidence limits as weighted means of the confidence
#' limits for the individual fits.
#'
#' @param object an object of class 'drc'.
#' @param fctList a list of non-linear functions to be compared.
#' @param respLev a numeric vector containing the response levels.
#' @param interval character string specifying the type of confidence intervals to be supplied.
#'   The default is "none". The choices "buckland" and "kang" are explained in the Details section.
#' @param linreg logical indicating whether or not additionally a simple linear regression model
#'   should be fitted.
#' @param clevel character string specifying the curve id in case on estimates for a specific
#'   curve or compound is requested. By default estimates are shown for all curves.
#' @param level numeric. The level for the confidence intervals. The default is 0.95.
#' @param type character string. Whether the specified response levels are absolute or
#'   relative (default).
#' @param display logical. If TRUE results are displayed. Otherwise they are not (useful in
#'   simulations).
#' @param na.rm logical indicating whether or not NA occurring during model fitting should be
#'   left out of subsequent calculations.
#' @param extended logical specifying whether or not an extended output (including fit summaries)
#'   should be returned.
#'
#' @return A matrix with two or more columns, containing the estimates and the corresponding
#'   estimated standard errors and possibly lower and upper confidence limits.
#'
#' @references Buckland, S. T. and Burnham, K. P. and Augustin, N. H. (1997) Model Selection:
#'   An Integral Part of Inference, \emph{Biometrics} \bold{53}, 603--618.
#'
#'   Kang, Seung-Ho and Kodell, Ralph L. and Chen, James J. (2000) Incorporating Model
#'   Uncertainties along with Data Uncertainties in Microbial Risk Assessment,
#'   \emph{Regulatory Toxicology and Pharmacology} \bold{32}, 68--72.
#'
#' @author Christian Ritz
#'
#' @seealso The function \code{\link{mselect}} provides a summary of fit statistics for
#'   several models fitted to the same data.
#'
#' @examples
#' ## Fitting an example dose-response model
#' ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
#'
#' ## Doing the actual model-averaging
#' maED(ryegrass.m1,
#' list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5)),
#' c(10, 50, 90))
#'
#' @keywords models nonlinear
"maED" <- function(object, fctList = NULL, respLev, interval = c("none", "buckland", "kang"), linreg = FALSE,
clevel = NULL, level = 0.95, type = c("relative", "absolute"), display = TRUE, na.rm = FALSE, extended = FALSE)
{
    interval <- match.arg(interval)
    type <- match.arg(type)
#    print(linreg)
#    print(level)
    
    ## Handling multiple curves in a single dataset
    ncolPM <- ncol(object$"parmMat")
    if ((!identical(ncolPM, 1)) && (is.null(clevel)))  # is also TRUE for a single curve!!!
    {
        retMat <- NULL
        for (i in 1:ncolPM)
        {
#            print((colnames(object$"parmMat"))[i])
            curveId <- (colnames(object$"parmMat"))[i]
#            cat(curveId, "\n") 
            retMat <- rbind(retMat, 
            maED(object, fctList, respLev, interval, linreg = linreg, clevel = curveId, level = level, 
            type = type, display = display, na.rm = na.rm, extended = extended))
        }
        return(retMat)
    } else {     # May 6 2010

    interval <- match.arg(interval)

    msMat <- do.call("mselect", list(object = object, fctList = fctList, sorted = "no")) 
       
#    expVec <- as.vector(exp(-msMat[, 2] / 2))
#    wVec <- expVec / sum(expVec, na.rm = na.rm)  
    # maybe better "combined" na.rm approach for edEst and wVec
    
#    ## Removing poor fits completely via a threshold (good approach?)
#    wVec[wVec < 0.01] <- 0
    
    lenfl <- length(fctList)
    lenrl <- length(respLev)
#    uniCID <- unique(S.alba.m1$data[, 4])
#    lenuniCID <- length(uniCID)
    
    numRows <- lenfl + 1
    numCols <- lenrl    
#    numCols <- lenrl * lenuniCID
        
    edEst <- matrix(NA, numRows + linreg, numCols)
    edSe <- matrix(NA, numRows + linreg, numCols)    
#    print(c(numRows, numRows + linreg))
#    print(linreg)
    
    ## Defining 'interval' argument for ED
    if (identical(interval, "kang"))
    {
        interval2 <- "delta"
    } else {
        interval2 <- "none"    
    }
    
    ## Calculating estimated ED values
#    edMat <- ED(object, respLev, interval2, clevel, type = type, display = FALSE, multcomp = TRUE)[["EDdisplay"]]
    edMat <- ED(object, respLev, interval2, clevel, type = type, display = FALSE)
    edEst[1, ] <- as.vector((edMat)[, 1])
    edSe[1, ] <- as.vector((edMat)[, 2])
    
    if (identical(interval2, "delta"))
    {
        edCll <- matrix(NA, numRows, numCols)
        edClu <- matrix(NA, numRows, numCols)
        
        edCll[1, ] <- as.vector((edMat)[, 3])
        edClu[1, ] <- as.vector((edMat)[, 4])        
    }
    for (i in 1:lenfl)
    {
#        edMati <- try(ED(update(object, fct = fctList[[i]]), respLev, interval2, clevel, type = type, display = FALSE, multcomp = TRUE)[["EDdisplay"]], silent = TRUE)
      edMati <- try(ED(update(object, fct = fctList[[i]]), respLev, interval2, clevel, 
                       type = type, display = FALSE), silent = TRUE)      
        if (inherits(edMati, "try-error"))
        {
            edMati <- matrix(NA, length(respLev), 4)
        }
         
        edEst[i + 1, ] <- as.vector((edMati)[, 1])
        edSe[i + 1, ] <- as.vector((edMati)[, 2])
        if (identical(interval2, "delta"))
        {
            edCll[i + 1, ] <- as.vector((edMati)[, 3])
            edClu[i + 1, ] <- as.vector((edMati)[, 4])
        }
    }
    
#    print(edEst)
    ## Adding simple linear regression fit
    if (linreg)
    {
        linFit1 <- lm(object$"data"[, 2:1])
        edLin <- ED.lin(linFit1, respLev)
        edEst[lenfl + 2, ] <- unlist((edLin)[, 1])
        edSe[lenfl + 2, ] <- unlist((edLin)[, 2])

        ## Updating weights
        expVec <- as.vector(exp(-c(msMat[, 2], AIC(linFit1)) / 2))
    
    } else {
        expVec <- as.vector(exp(-msMat[, 2] / 2))
    }
#    print(edEst)
    wVec <- expVec / sum(expVec, na.rm = na.rm)  

    edVec <- apply(edEst * wVec, 2, sum, na.rm = na.rm)
    if (identical(interval, "none"))
    {
        retMat <- as.matrix(cbind(edVec))
        colnames(retMat) <- colnames(edMat)[1]
    }
    if (identical(interval, "buckland"))
    {
        seVec <- apply(sqrt(edSe^2 + (t(t(edEst) - apply(edEst, 2, mean, na.rm = na.rm)))^2) * wVec, 2, 
        sum, na.rm = na.rm)
### Thresholding                
#        iVec <- wVec < 0.01        
#        seVec <- apply(sqrt(edSe[iVec, ]^2 + (t(t(edEst[iVec, ]) - apply(edEst[iVec, ], 2, 
#        mean, na.rm = na.rm)))^2) * wVec[iVec], 2, sum, na.rm = na.rm)
        quantVal <- qnorm(1 - (1 - level)/2) * seVec
        retMat <- as.matrix(cbind(edVec, seVec, edVec - quantVal, edVec + quantVal))
        colnames(retMat) <- c(colnames(edMat)[c(1, 2)], "Lower", "Upper")
    }
    if (identical(interval, "kang"))
    {  
        retMat <- as.matrix(cbind(apply(edEst * wVec, 2, sum, na.rm = na.rm), 
        apply(edCll * wVec, 2, sum, na.rm = na.rm), 
        apply(edClu * wVec, 2, sum, na.rm = na.rm)))
        colnames(retMat) <- colnames(edMat)[c(1,3,4)]
    }   
    rownames(retMat) <- rownames(edMat)

    ## Constructing matrix of fit summaries 
    disMat <- as.matrix(cbind(edEst, wVec))
#    colnames(disMat) <- c(paste("EC", rownames(edMat), sep = ""), "Weight")
    colnames(disMat) <- c(paste("ED", respLev, sep = ""), "Weight")
#    rownames(disMat) <- rownames(msMat)
    if (linreg) 
    {
        rownames(disMat) <- c(rownames(msMat), "Lin")
    } else {
        rownames(disMat) <- rownames(msMat)
    }

#    if (lenuniCID > 1)
#    {
#        rownames(disMat) <- paste(rownames(msMat), uniCID, sep = "-")
#    } else {
#        rownames(disMat) <- rownames(msMat)    
#    }
    if (display)
    {
        print(disMat)
        cat("\n") 
    }
    
#    resPrint(resMat, "Estimated effective doses", interval, "Model-averaging", display)
    if (extended)
    {
        return(list(estimates = retMat, fits = disMat))
    } else {
        retMat
    }
    }   #May 6 2010
}