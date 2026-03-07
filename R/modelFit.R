#' @title Assessing the model fit
#'
#' @description
#' Checking the fit of a dose-response model by means of formal significance tests.
#'
#' @param object object of class 'drc'.
#' @param test character string defining the test method to apply.
#' @param method character string specifying the method to be used for assessing the model fit.
#'
#' @details
#' Currently two methods are available. For continuous data the classical lack-of-fit test is
#' applied (Bates and Watts, 1988). The test compares the dose-response model to a more general
#' ANOVA model using an approximate F-test. For quantal data the crude goodness-of-fit test
#' based on Pearson's statistic is used.
#'
#' None of these tests are very powerful. A significant test result is more alarming than a
#' non-significant one.
#'
#' @return An object of class 'anova' which will be displayed in much the same way as an
#'   ordinary ANOVA table.
#'
#' @references
#' Bates, D. M. and Watts, D. G. (1988)
#' \emph{Nonlinear Regression Analysis and Its Applications},
#' New York: Wiley & Sons (pp. 103--104).
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Comparing the four-parameter log-logistic model
#' ##  to a one-way ANOVA model using an approximate F test
#' ## in other words applying a lack-of-fit test
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
#' modelFit(ryegrass.m1)
#'
#' @keywords models nonlinear
"modelFit" <- function(object, test = NULL, method = c("gof", "cum"))
{
    method <- match.arg(method)

    ## Fitting ANOVA model
    testList <- switch(object$"type", 
    "continuous" = drmLOFls(), 
    "binomial" = drmLOFbinomial(), 
    "Poisson" = drmLOFPoisson())
  
    switch(object$"type", 
    binomial = gofTest(object, testList$"gofTest"),  
    continuous = lofTest(object, testList$"anovaTest"))
}


"lofTest" <- function(object, anovaTest)
{
    if (!is.null(anovaTest))
    {
        dose <- object$"dataList"$"dose"
        resp <- object$"dataList"$"resp"
        curveid <- object$"dataList"$"curveid"
        if (is.null(object$"boxcox"))
        {
            bcAdd <- 0
        } else {
            bcAdd <- object$"boxcox"$"bcAdd"
        }
        afList <- anovaFormula(dose, resp, curveid, bcAdd)
        anovaForm <- afList$"anovaFormula"
        anovaData <- afList$"anovaData"        
        anovaModel <- anovaTest(anovaForm, anovaData)
        if (is.null(anovaModel))
        {
            return(returnFct())
        }

        anovaDF <- df.residual(anovaModel$"anovaFit")
        nlsDF <- df.residual(object)
        dfModel <- c(anovaDF, nlsDF)
        dfDiff <- c(NA, (nlsDF - anovaDF))
                                
        if (identical(anovaModel$"test", "F"))
        {
            anovaSS <- deviance(anovaModel$"anovaFit")
            anovaDF <- df.residual(anovaModel$"anovaFit")
            dfModel <- c(anovaDF, nlsDF)
            nlsSS <- object$"fit"$"value"
            loglik <- c(anovaSS, nlsSS)

            testStat <- (nlsSS - anovaSS)/dfDiff[2]/(anovaSS/anovaDF)
            pVal <- c(NA, pf(testStat, dfDiff[2], anovaDF, lower.tail = FALSE))
            testStat <- c(NA, testStat)

            headName<-"Lack-of-fit test\n"
            rowNames<-c("ANOVA", "DRC model")
            colNames<-c("ModelDf", "RSS", "Df", "F value", "p value")    
        }
            
        if (identical(anovaModel$"test", "lr"))
        {
            anovaDF <- anovaDF + (object$"sumList"$"lenData" - dim(anovaModel$"anovaFit"$"data")[1])
            dfModel <- c(anovaDF, nlsDF)
            dfDiff <- c(NA, (nlsDF - anovaDF))    
            loglik <- c(logLik(anovaModel$"anovaFit"), logLik(object))

            testStat <- 2*(loglik[1] - loglik[2])
            pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
            testStat <- c(NA, testStat)

            headName <- "Goodness-of-fit test\n"
            rowNames <- c("ANOVA", "DRC model")
            colNames <- c("ModelDf", "Log lik", "Df", "Chisq value", "p value")                
        }
        return(returnFct(dfModel, loglik, dfDiff, testStat, pVal, headName, colNames, rowNames))
         
    } else {
        return(returnFct())
    }
}

"gofTest" <- function(object, gofTest)
{       
    gofTest <- gofTest(object$"dataList"$"resp", weights(object), fitted(object), df.residual(object))
    
    if (!is.null(gofTest))
    {
        returnFct(c(NA, NA), c(NA, NA), c(NA, gofTest[2]), c(NA, gofTest[1]), 
        c(NA, 1 - pchisq(gofTest[1], gofTest[2])), "Goodness-of-fit test\n",
        c("", "", "Df", "Chisq value", "p value"), c("", "DRC model"))
    } else {  
        returnFct()
    }
}


"returnFct" <- function(dfModel = c(NA, NA), loglik = c(NA, NA), dfDiff = c(NA, NA), testStat = c(NA, NA), 
pVal = c(NA, NA), headName = "No test available\n", colNames = c("ModelDf", "Log lik", "Df", "Chisq value", "p value"), 
rowNames = c("", "DRC model"))
{
    dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
    dimnames(dataFra) <- list(rowNames, colNames)
    structure(dataFra, heading = headName, class = c("anova", "data.frame"))  
}