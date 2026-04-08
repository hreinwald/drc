#' Simulation of dose-response data and ED estimation
#'
#' Simulates dose-response datasets using parametric or non-parametric methods and estimates
#' effective doses (ED values) from each simulated dataset. Useful for assessing the
#' performance of ED estimation methods via Monte Carlo simulation.
#'
#' @param noSim integer. Number of simulations to run.
#' @param edVal numeric vector of ED levels to estimate (default is \code{c(10, 20, 50)}).
#' @param type character string. Either "non-parametric" or "parametric" simulation.
#' @param response character string. Either "bin" (binomial) or "con" (continuous) response.
#' @param fct dose-response function used for simulation (default is \code{LL.2()}).
#' @param coefVec numeric vector of model coefficients for parametric simulation.
#' @param method character string. Estimation method: "sp" (semi-parametric), "p" (parametric),
#'   or "np" (non-parametric).
#' @param doseVec numeric vector of dose values.
#' @param nVec numeric vector of sample sizes per dose (for binomial response).
#' @param pVec numeric vector of expected response probabilities (for non-parametric simulation).
#' @param rVec numeric vector of responses.
#' @param resVar numeric. Residual variance (for continuous response).
#' @param pfct dose-response function used for fitting (defaults to \code{fct}).
#' @param reference character string specifying the reference for ED estimation.
#' @param span numeric. Smoothing parameter for local regression. NA uses default.
#' @param minmax character string. Type of min/max calculation. Default is "response".
#' @param lower numeric. Lower bounds for optimization.
#' @param upper numeric. Upper bounds for optimization.
#' @param seedVal integer. Random seed for reproducibility (default is 200810201).
#'
#' @return A list with components \code{edArray} (array of ED estimates), \code{mixVec},
#'   \code{edVal}, \code{aicVec}, and \code{spanVec}.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"simFct" <- function(noSim, edVal = c(10, 20, 50), type = c("non-parametric", "parametric"), 
response = c("bin", "con"), fct = LL.2(), coefVec, method = c("sp", "p", "np"), 
doseVec, nVec, pVec, rVec, resVar, pfct = fct, reference = NULL, span = NA, 
minmax = "response", lower = NULL, upper = NULL, seedVal = 200810201)
{
    method <- match.arg(method)
    response <- match.arg(response)
    type <- match.arg(type)

    set.seed(seedVal)

    lenData <- length(doseVec)  # replace lenpv?

    ## Parametric simulations
    ## Drawing random dose-response curves
    if (type == "parametric")
    {
        ## Model fit to simulate from
        if (response == "bin")
        {
            simMat <- rdrm(noSim, fct, coefVec, doseVec, yerror = "rbinom", ypar = nVec, onlyY = TRUE)
        } else {
            simMat <- rdrm(noSim, fct, coefVec, doseVec, ypar = c(0, sqrt(resVar)), onlyY = TRUE)
        }       
        ## drop = FALSE also in rdrm???
        print(simMat$y[1, ])        
    }

    ## Non-parametric simulations
    if (type == "non-parametric")
    {
        lenpv <- length(pVec)
    
        simMat <- matrix(NA, noSim, lenpv)
        if (response == "bin")
        {
            for (i in 1:noSim)
            {
                simMat[i, ] <- rbinom(lenpv, nVec, pVec)
            }
        } else {
            for (i in 1:noSim)
            {
                simMat[i, ] <- rnorm(lenpv, pVec, sqrt(resVar))
            }
        }
        simMat <- list(y = simMat)
        print(simMat$y[1, ])
    }
    

    lenev <- length(edVal)
    aicVec <- rep(NA, noSim)
    edMat <- array(NA, c(lenev, 3, noSim))
    mixVec <- rep(NA, noSim)
    spanVec <- rep(span, noSim)

    ## Generalized cross-validation criterion
    gcvFct <- function(doseVec, y)  # define function outside the i loop
    {
        gcvVec <- rep(NA, 20)
        for (j in 1:20)
        {
            tempLoess <- try(loess(y ~ doseVec, degree = 1, span = j/20), silent = TRUE)
            if (inherits(tempLoess, "try-error"))
            {
                gcvVec[j] <- NA
            } else {
                gcvVec[j] <- sum(residuals(tempLoess)^2) / (lenData - tempLoess$trace.hat)^2
            }
        }
        ((1:20)/20)[which.min(gcvVec)]
    }

    for (i in 1:noSim) 
    {
        ## Converting to proportions
        if (response == "bin")
        {
            y <- simMat$y[i, ] / nVec
        } else {
            y <- simMat$y[i, ]
        }
        
        ## Obtaining model-robust fit 
        if (method == "sp")
        {      
            parModel <- drm(y ~ doseVec, fct = pfct)
            
            if (is.na(span)) 
            {
                spanVec[i] <- gcvFct(doseVec, y)
            }
            loessModel <- loess(y ~ doseVec, degree = 1, span = spanVec[i])

            tempModel <- mrdrm(parModel, loessModel)
            if (inherits(tempModel, "try-error"))
            {
                tempModel <- list(edMat = NA, mixing = NA, aic = NA)
            }  else {
                 aicVec[i] <- tempModel$gof[3]
                 mixVec[i] <- tempModel$lambda                 
                 edMat[, , i] <- ED(tempModel, edVal, interval = "approximate", minmax = minmax, 
                 lower = lower, upper = upper, display = FALSE)[, c(1:3)]
            }
        }
        if (method == "p")
        {
            if (response == "con")
            {
                tempModel <- try(drm(formula = y ~ doseVec, fct = pfct), silent = TRUE)
            } else {
                tempModel <- try(drm(formula = y ~ doseVec, weights = nVec, fct = pfct, type = "binomial"), 
                silent = TRUE)
            }
            if (inherits(tempModel, "try-error"))
            {
                edMat[, , i] <- NA 
                mixVec[i] <- NA
            } else {
                tempED <- try(ED(tempModel, edVal, display = FALSE, interval = "delta")[, c(1, 3, 4)], silent = TRUE)
                if (inherits(tempED, "try-error"))
                {
                    edMat[, , i] <- NA
                    mixVec[i] <- NA
                } else {
                    edMat[, , i] <- tempED
                    mixVec[i] <- 0
                    aicVec[i] <- AIC(tempModel)
                }
            }
        }     
    }
    list(edArray = edMat, mixVec = mixVec, edVal = edVal, aicVec = aicVec, spanVec = spanVec)
}

## Calculating coverage percentage
coverFct <- function(mfit, simres, edVec = NULL)
{
    edVal <- simres$edVal
    if (is.null(edVec)) 
    {
        edVec <- ED(mfit, edVal, display = FALSE)[, 1]
    }

    lenem <- length(edVal)
    cpVec <- rep(NA, lenem)
    cplVec <- rep(NA, lenem)
    cpuVec <- rep(NA, lenem)        
    mvVec <- rep(NA, lenem)    
    mwVec <- rep(NA, lenem)
    notNA <- rep(NA, lenem)    
    names(cpVec) <- edVal
    for (i in 1:lenem)
    {
        notNA[i] <- sum((!is.na(simres$edArray[i, 2, ])) & (!is.na(simres$edArray[i, 3, ])))
        cplVec[i] <- sum(is.na(simres$edArray[i, 2, ]) & (simres$edArray[i, 3, ] > edVec[i]), na.rm = TRUE)
        cpuVec[i] <- sum(is.na(simres$edArray[i, 3, ]) & (simres$edArray[i, 2, ] < edVec[i]), na.rm = TRUE)        
        cpVec[i] <- sum((simres$edArray[i, 2, ] < edVec[i]) & (simres$edArray[i, 3, ] > edVec[i]), na.rm = TRUE) / notNA[i]
        mvVec[i] <- mean(simres$edArray[i, 1, ], na.rm = TRUE)
        mwVec[i] <- mean(simres$edArray[i, 3, ] - simres$edArray[i, 2, ], na.rm = TRUE)
    }
    list(coverage = cpVec, covLow = cplVec, covUp = cpuVec, true = edVec, mean = mvVec, width = mwVec, 
    notNAs = notNA, NAs = length(simres$edArray[1, 2,]) - notNA, mixingAverage = mean(simres$mixVec, na.rm = TRUE))
}
