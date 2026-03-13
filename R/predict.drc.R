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
#        predValues <- fitted(object)  # not used
#        newdata <- data.frame(object$data[, 1], object$data[, 3])
#        dataList <- object[["dataList"]]

        ## New part (25/6-2014)
        doseVec <- dataList[["dose"]]
        if (identical(respType, "event"))
        {
            groupLevels <- as.character(dataList[["plotid"]])
        } else {
            groupLevels <- as.character(dataList[["curveid"]])
        }        
#                
#        if (identical(respType, "event"))
#        {
#            newdata <- data.frame(dataList[["dose"]], dataList[["plotid"]])
#        } else {
#            newdata <- data.frame(dataList[["dose"]], dataList[["curveid"]])
#        }
    } else {
        
        if (checkND)
        {  
            dName <- dataList[["names"]][["dName"]]
            if (any(names(newdata) %in% dName))
            {
                doseVec <- newdata[, dName]  
            } else {
                doseVec <- newdata[, 1]
#                warning("Dose variable not in 'newdata'")
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
            ndRows <- if (is.data.frame(newdata) || is.matrix(newdata)) nrow(newdata) else length(newdata)
            groupLevels <- rep(1, ndRows)
        }
#         
#         
#         if (ncol(newdata) < (doseDim + 1)) {newdata <- data.frame(newdata, rep(1, nrow(newdata)))}
# #        ndncol <- ncol(newdata)
# #        doseVec <- newdata[, 1:(ndncol-1)]
#         doseVec <- newdata[, 1:doseDim]
# #        groupLevels <- as.character(newdata[, ndncol])  # 'as.character()' used to suppress factor levels         
#         groupLevels <- as.character(newdata[, doseDim + 1])  # 'as.character()' used to suppress factor levels         
    }
    noNewData <- length(groupLevels)
    
#    if (ncol(newdata) < 2) {newdata <- data.frame(newdata, rep(1, nrow(newdata)))}
#    if (ncol(newdata) > 2) {stop("More than 2 variables in 'newdata' argument")}
    
    ## Defining dose values -- dose in the first column!   
#    doseVec <- newdata[, 1]
#    groupLevels <- as.character(newdata[, 2])  # 'as.character()' used to suppress factor levels 
#    noNewData <- length(doseVec)


    ## Transforming to dose scale if necessary
    powerExp <- (object$"curve")[[2]]
    if (!is.null(powerExp))
    {
        doseVec <- powerExp ^ doseVec
    }

    ## Retrieving matrix of parameter estimates
    parmMat <- object[["parmMat"]] 
    pm <- t(parmMat[, groupLevels, drop = FALSE])
           
#    parmNames <- colnames(parmMat)
#    lenCN <- length(parmNames)
#    indVec <- 1:lenCN
#    names(indVec) <- parmNames
#    if (lenCN > 1)
#    {
#        indVec <- indVec[as.character(newdata[, 2])]
#        
##        groupLevels <- newdata[, 2]
#        if (!all(is.numeric(groupLevels)))
#        {
##            pm <- parmMat[, as.character(groupLevels)]  # 'as.character()' used to suppress factor levels            
#            pm <- parmMat[, groupLevels]            
#        } else {
#            pm <- parmMat[, groupLevels]
#        }
#        pm <- parmMat[, groupLevels]
#        
#    } else {
#        lenDV <- length(doseVec)
##        indVec <- rep(1, lenDV)
#        pm <- matrix(parmMat[, 1], length(parmMat[, 1]), lenDV)
#    }    

#    ## Checking for NAs in matrix of parameter estimates
#    naVec <- rep(FALSE, lenCN)
#    for (i in 1:lenCN)
#    {
#        naVec[i] <- any(is.na(parmMat[, i]))
#    }
#    parmMat <- parmMat[, !naVec, drop = FALSE] 


    ## Retrieving variance-covariance matrix
    sumObj <- summary(object, od = od)
#    varMat <- sumObj[["varMat"]]  
    vcovMat <- vcov.(object)      

    ## Defining index matrix for parameter estimates
    indexMat <- object[["indexMat"]]
    
    ## Calculating predicted values  
#    indexVec <- as.vector(indVec)  
#    print(indexVec)  
#    lenIV <- length(indexVec)    
    
    
#    retMat <- matrix(0, lenIV, 4)
    retMat <- matrix(0, noNewData, 4)
    colnames(retMat) <- c("Prediction", "SE", "Lower", "Upper")
    objFct <- object[["fct"]]
#    print(pm)
#    print(doseVec)
    retMat[, 1] <- objFct$"fct"(doseVec, pm)
#    print(pm)
    
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
#                lmObj <- lm(seVec ~ estVec)  # linear, not great
                lmObj <- lm(log(seVec) ~ log(estVec))
                sePred <- exp(predict(lmObj, data.frame(estVec = doseVec)))
            } else {
                sePred <- ssdSEfct(estVec, seVec, doseVec)
            }
#            print(sePred)
#            print(object[["fct"]][["derivx"]](doseVec, pm))
            derivxRes <- object[["fct"]][["derivx"]](doseVec, pm)
#            print(derivxRes)
            # if (is.finite(derivxRes))
            # {
            #     sumObjRV <- (derivxRes * sePred)^2  
            # } else {
            #     sumObjRV <- 0  # setting Inf * 0 = 0 (would be NaN otherwise) 
            # } 
            sumObjRV <- rep(0, length(derivxRes))
            isFinDR <- is.finite(derivxRes) 
            sumObjRV[isFinDR] <- ((derivxRes * sePred)^2)[isFinDR]
            
#            print(sumObjRV)
        } 
        if (identical(interval, "prediction"))
        {
            sumObjRV <- rep(sumObj$"resVar", noNewData)
        } 
        #else {
        #     sumObjRV <- 0
        #  }
#        rowIndex <- 1    
#        for (i in indexVec)
#        for (i in 1:ncol(indexMat))

#        groupLevels <- newdata[, 2]
        piMat <- indexMat[, groupLevels, drop = FALSE]
#        print(piMat)
#        print(groupLevels)
        for (rowIndex in 1:noNewData)
        {
#            parmInd <- indexMat[, i]
#            print(indexVec) 
#            print(varMat)
#            print(parmInd)       
            
#            varCov <- varMat[parmInd, parmInd]
#            print(varCov)
#            groupLevels <- newdata[, 2]
#            parmInd <- indexMat[, groupLevels[rowIndex]]
#            varCov <- varMat[parmInd, parmInd]

            parmInd <- piMat[, rowIndex] 
            varCov <- vcovMat[parmInd, parmInd]

#            parmChosen <- t(parmMat[, i, drop = FALSE])
#            parmChosen <- t(pm[, rowIndex, drop = FALSE])
#            dfEval <- deriv1(doseVec[rowIndex], parmChosen)

            dfEval <- deriv1(doseVec[rowIndex], pm[rowIndex, , drop = FALSE])
            varVal <- dfEval %*% varCov %*% dfEval
            retMat[rowIndex, 2] <- sqrt(varVal)  
#            retMat[rowIndex, 2] <- sqrt(dfEval %*% varCov %*% dfEval)  

            if (!se.fit)
            {
                #retMat[rowIndex, 3:4] <- rep(retMat[rowIndex, 1], 2) + 
                #  (tquan * sqrt(varVal + sumObjRV[rowIndex])) * c(-1, 1)
                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObjRV[rowIndex])
                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObjRV[rowIndex])   
            }    
#            if (identical(interval, "confidence"))
#            {
#                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal)
#                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal)            
#            }
#            if (identical(interval, "prediction"))
#            {
#                sumObjRV <- sumObj$"resVar"
#                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObjRV)
#                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObjRV)                        
#            }          
#            rowIndex <- rowIndex + 1        
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


