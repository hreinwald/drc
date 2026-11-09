#' @title Interactive dose-response modelling
#'
#' @description Internal helper used by \code{\link{drm}} when \code{separate = TRUE}: the model is
#' fitted to each curve separately and the fits are combined into one \code{drc} object.
#'
#' @param x dose values.
#' @param y response values.
#' @param curveid curve identifiers.
#' @param weights weights.
#' @param fct model function (list).
#' @param type response type.
#' @param control list from \code{\link{drmc}}.
#' @param robust,bcVal,bcAdd,lowerl,upperl,logDose,start arguments of \code{\link{drm}} that are
#'   forwarded to the per-curve fits.
#'
#' @keywords internal
"idrm" <- function(x, y, curveid, weights, fct, type, control, robust = "mean", bcVal = NULL, 
                   bcAdd = 0, lowerl = NULL, upperl = NULL, logDose = NULL, start = NULL)
{
    oneFunction <- !is.list(fct[[1]])
 
    ## Fitting models for each curve
    fitList <- list() 
    uniCur <- unique(curveid)   
    numCur <- length(uniCur)
    for (i in 1:numCur)
    {
        if (oneFunction)
        {
            if (is.null(start))
            {
                fitList[[i]] <- drm(y~x, curveid, weights = weights, fct = fct, type = type, 
                subset = curveid == uniCur[i], separate = FALSE, control = control,
                robust = robust, bcVal = bcVal, bcAdd = bcAdd, lowerl = lowerl, upperl = upperl, 
                logDose = logDose)
            } else {
                fitList[[i]] <- drm(y~x, curveid, weights = weights, fct = fct, type = type, 
                subset = curveid == uniCur[i], separate = FALSE, control = control,
                robust = robust, bcVal = bcVal, bcAdd = bcAdd, lowerl = lowerl, upperl = upperl, 
                logDose = logDose, start = start)
            }
        } else {
            tempFitlist <- list()
            for (j in 1:length(fct))
            {
                tempFitlist[[j]] <- drm(y~x, curveid, weights = weights, fct = fct[[j]], 
                type = type, subset = curveid == uniCur[i], separate = FALSE, control = control)
            }
            fitList[[i]] <- tempFitlist
        }
    }

    retList <- fitList[[1]]
    if (oneFunction)
    {
        dataMat <- fitList[[1]]$"data"   
        tdataList <- fitList[[1]]$"dataList"   
        
        parmMat <- fitList[[1]]$"parmMat"
        nlsFit <- fitList[[1]]$"nlsFit"
        
        parNames <- fitList[[1]]$"parNames"
        numPar <- length(parNames[[1]])        
        parNames[[3]] <- rep(uniCur[i], numPar) 
               
        if (numCur > 1)
        {
            nlsFit[[1]] <- nlsFit
            pnList <- list()
            pnList[[1]] <- parNames
            pnList[[1]][[3]] <- rep(uniCur[1], numPar)
            
            for (i in 2:numCur)
            {
                parmMat <- cbind(parmMat, fitList[[i]]$"parmMat")
                dataMat <- rbind(dataMat, fitList[[i]]$"data")  
                tdataList <- mapply(c, tdataList, fitList[[i]]$"dataList", SIMPLIFY = FALSE)              
                nlsFit[[i]] <- fitList[[i]]$"nlsFit"
                pnList[[i]] <- fitList[[i]]$"parNames"
                pnList[[i]][[3]] <- rep(uniCur[i], numPar)
            }
        }
        retList$"dataList" <- tdataList 
        retList$"dataList"$"names" <- fitList[[1]]$"dataList"$"names"        
        retList$"data" <- dataMat 
        retList$"parmMat" <- parmMat 
        
        plotFct <- function(x) {matrix(unlist(lapply(fitList, function(y)y$"curve"[[1]](x))), ncol = numCur)}
        retList$"curve" <- list(plotFct, fitList[[1]]$"curve"[[2]])

        bVec <- as.vector(unlist(lapply(pnList, function(x){x[[2]]})))
        cVec <- as.vector(unlist(lapply(pnList, function(x){x[[3]]}))) 
        aVec <- paste(bVec, cVec, sep = ":")
        retList$"parNames" <- list(aVec, bVec, cVec)       
    
        retList$"indexMat" <- matrix(c(1:(numCur * numPar)), numPar, numCur)
        names(fitList) <- uniCur
        retList$"objList" <- fitList
    
        coefVec <- as.vector(unlist(lapply(fitList, function(x){x$"fit"$"par"})))
        names(coefVec) <- aVec
        retList$"coefficients" <- coefVec
        
        retList$"df.residual" <- sum(unlist(lapply(fitList, function(x){x$"df.residual"})))
        retList$"minval" <- sum(unlist(lapply(fitList, function(x){x$"fit"$"value"})))
        
        retList$"fit" <- nlsFit
    
        
    } else {
    
    }
        
    class(retList) <- c("drc")
    return(retList)
}


