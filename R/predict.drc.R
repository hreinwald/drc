#' @title Prediction
#'
#' @description
#' Predicted values for models of class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param newdata an optional data frame in which to look for variables with
#'   which to predict. If omitted, the fitted values are used.
#' @param se.fit logical. If TRUE standard errors are required.
#' @param interval character string. Type of interval calculation:
#'   \code{"none"}, \code{"confidence"}, \code{"prediction"}, or \code{"ssd"}.
#' @param level tolerance/confidence level.
#' @param na.action function determining what should be done with missing values
#'   in \code{newdata}. The default is to predict \code{NA}.
#' @param od logical. If TRUE adjustment for over-dispersion is used.
#' @param vcov. function providing the variance-covariance matrix.
#'   \code{\link{vcov}} is the default, but \code{sandwich} is also an option
#'   (for obtaining robust standard errors).
#' @param ssdSEfct specifies the function for interpolating standard errors
#'   between observed standard errors. The default is linear interpolation on
#'   log-log scale (back-transformed).
#' @param constrain logical. If TRUE (default) predicted values are truncated
#'   within meaningful limits, i.e., 0 and, possibly, 1.
#' @param checkND logical indicating whether or not names in \code{newdata}
#'   data frame match the names in the original data frame used for fitting
#'   the model. Default is TRUE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A matrix with as many rows as there are dose values provided in
#'   \code{newdata} or in the original dataset (in case \code{newdata} is not
#'   specified) and, at most, 4 columns containing fitted values, standard
#'   errors, lower and upper limits of confidence/prediction intervals.
#'
#' @seealso For details see the help page for \code{\link{predict.lm}}.
#'
#' @examples
#' ## Fitting a model
#' spinach.model1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
#'
#' ## Predicting values at dose=2 (with standard errors)
#' predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")), se.fit = TRUE)
#'
#' ## Getting confidence intervals
#' predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")),
#' interval = "confidence")
#'
#' ## Getting prediction intervals
#' predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")),
#' interval = "prediction")
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"predict.drc" <- function(object, newdata, se.fit = FALSE, 
                          interval = c("none", "confidence", "prediction", "ssd"), 
                          level = 0.95, na.action = na.pass, od = FALSE, vcov. = vcov, 
                          ssdSEfct = NULL, constrain = TRUE, checkND = TRUE, ...)
{
    ## Checking arguments
    interval <- match.arg(interval)
    respType <- object[["type"]]

    dataList <- object[["dataList"]]    
    doseDim <- ncol(dataList[["dose"]])
    if (is.null(doseDim)) {doseDim <- 1} 

    ## Assigning dataset from object if no data frame is provided
    if (missing(newdata)) 
    {
        ## New part (25/6-2014)
        doseVec <- dataList[["dose"]]
        if (identical(respType, "event"))
        {
            groupLevels <- as.character(dataList[["plotid"]])
        } else {
            groupLevels <- as.character(dataList[["curveid"]])
        }        
    } else {
        
        if (checkND)
        {  
            dName <- dataList[["names"]][["dName"]]
            if (any(names(newdata) %in% dName))
            {
                doseVec <- newdata[, dName]  
            } else {
                doseVec <- newdata[, 1]
            }
        } else {
            doseVec <- newdata
        }
        
        cName <- dataList[["names"]][["cNames"]]
        if (any(names(newdata) %in% cName))
        {
            groupLevels <- as.character(newdata[, cName])  
            # as.character() removes factor encoding  
          
        } else {
            nRows <- if (is.data.frame(newdata) || is.matrix(newdata)) nrow(newdata) else length(newdata)
            groupLevels <- rep(1, nRows)
        }
    }
    noNewData <- length(groupLevels)
    
    ## Transforming to dose scale if necessary
    powerExp <- (object$"curve")[[2]]
    if (!is.null(powerExp))
    {
        doseVec <- powerExp ^ doseVec
    }

    ## Retrieving matrix of parameter estimates
    parmMat <- object[["parmMat"]] 
    pm <- t(parmMat[, groupLevels, drop = FALSE])

    ## Retrieving variance-covariance matrix
    sumObj <- summary(object, od = od)
    vcovMat <- vcov.(object)      

    ## Defining index matrix for parameter estimates
    indexMat <- object[["indexMat"]]
    
    ## Calculating predicted values  
    retMat <- matrix(0, noNewData, 4)
    colnames(retMat) <- c("Prediction", "SE", "Lower", "Upper")
    objFct <- object[["fct"]]
    retMat[, 1] <- objFct$"fct"(doseVec, pm)
    
    ## Checking if derivatives are available
    deriv1 <- objFct$"deriv1"
    if (is.null(deriv1))
    {
        return(retMat[, 1])        
    }    

    ## Calculating the quantile to be used in the confidence intervals
    if (!identical(interval, "none"))
    {    
        if (identical(respType, "continuous"))
        {
            tquan <- qt(1 - (1 - level)/2, df.residual(object))   
        } else {
            tquan <- qnorm(1 - (1 - level)/2)
        }
    }  
    
    ## Calculating standard errors and/or confidence intervals
    if (se.fit || (!identical(interval, "none")))
    {
        sumObjRV <- rep(0, noNewData)
        if (identical(interval, "ssd") & identical(object[["type"]], "ssd"))
        {
            estVec <- object[["dataList"]][["dose"]]
            seVec <- object[["dataList"]][["weights"]]
            
            if (is.null(ssdSEfct)) 
            {
                lmObj <- lm(log(seVec) ~ log(estVec))
                sePred <- exp(predict(lmObj, data.frame(estVec = doseVec)))
            } else {
                sePred <- ssdSEfct(estVec, seVec, doseVec)
            }
            derivxRes <- object[["fct"]][["derivx"]](doseVec, pm)
            # if (is.finite(derivxRes))
            # {
            #     sumObjRV <- (derivxRes * sePred)^2  
            # } else {
            #     sumObjRV <- 0  # setting Inf * 0 = 0 (would be NaN otherwise) 
            # } 
            sumObjRV <- rep(0, length(derivxRes))
            isFinDR <- is.finite(derivxRes) 
            sumObjRV[isFinDR] <- ((derivxRes * sePred)^2)[isFinDR]
        } 
        if (identical(interval, "prediction"))
        {
            sumObjRV <- rep(sumObj$"resVar", noNewData)
        } 
        piMat <- indexMat[, groupLevels, drop = FALSE]
        for (rowIndex in 1:noNewData)
        {
            parmInd <- piMat[, rowIndex] 
            varCov <- vcovMat[parmInd, parmInd]

            dfEval <- deriv1(doseVec[rowIndex], pm[rowIndex, , drop = FALSE])
            varVal <- dfEval %*% varCov %*% dfEval
            retMat[rowIndex, 2] <- sqrt(varVal)  

            if (!se.fit)
            {
                #retMat[rowIndex, 3:4] <- rep(retMat[rowIndex, 1], 2) + 
                #  (tquan * sqrt(varVal + sumObjRV[rowIndex])) * c(-1, 1)
                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObjRV[rowIndex])
                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObjRV[rowIndex])   
            }    
        }
    }
    ## Imposing constraints on predicted values
    if (constrain)
    {
        objType <- object[["type"]]
        if (identical(objType, "binomial") || identical(objType, "ssd"))
        {
            retMat[, 4] <- pmin(retMat[, 4], 1)
        }
        if (!identical(objType, "continuous"))
        {
            retMat[, 3] <- pmax(retMat[, 3], 0)
        }
    }
    
    ## Keeping relevant indices
    keepInd <- 1
    if (se.fit) {keepInd <- c(keepInd, 2)}
    if (!identical(interval, "none")) {keepInd <- c(keepInd, 3, 4)}
    
    if (length(keepInd) > 1) {
        return(retMat[, keepInd, drop = FALSE])
    } else {
        return(retMat[, keepInd])
    }
}


