#' @title Fitting dose-response models
#'
#' @description A general model fitting function for analysis of various types of dose-response data.
#'
#' @param formula a symbolic description of the model to be fit. Either of the form
#'   \code{response ~ dose} or as a data frame with response values in first column and dose
#'   values in second column.
#' @param curveid a numeric vector or factor containing the grouping of the data.
#' @param pmodels a data frame with as many columns as there are parameters in the non-linear
#'   function. Or a list containing a formula for each parameter in the nonlinear function.
#' @param weights a numeric vector containing weights. For continuous/quantitative responses,
#'   inverse weights are multiplied inside the squared errors (weights should have the same unit
#'   as the response). For binomial responses weights provide information about the total number
#'   of binary observations used to obtain the response.
#' @param data an optional data frame containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the
#'   fitting process.
#' @param fct a list with three or more elements specifying the non-linear function, the
#'   accompanying self starter function, the names of the parameters in the non-linear function
#'   and, optionally, the first and second derivatives as well as information used for
#'   calculation of ED values. Use \code{\link{getMeanFunctions}} for a full list.
#' @param type a character string specifying the distribution of the data. The default is
#'   \code{"continuous"}, corresponding to a normal distribution. Other choices include
#'   \code{"binomial"}, \code{"Poisson"}, \code{"negbin1"}, \code{"negbin2"}, \code{"event"},
#'   and \code{"ssd"}.
#' @param bcVal a numeric value specifying the lambda parameter to be used in the Box-Cox
#'   transformation.
#' @param bcAdd a numeric value specifying the constant to be added on both sides prior to
#'   Box-Cox transformation. The default is 0.
#' @param start an optional numeric vector containing starting values for all mean parameters
#'   in the model. Overrules any self starter function.
#' @param na.action a function for treating missing values (\code{NA}s). Default is
#'   \code{\link{na.omit}}.
#' @param robust a character string specifying the rho function for robust estimation.
#'   Default is non-robust least squares estimation (\code{"mean"}). Available robust methods
#'   are: \code{"median"}, \code{"lms"}, \code{"lts"}, \code{"trimmed"}, \code{"winsor"}, and
#'   \code{"tukey"}.
#' @param logDose a numeric value or \code{NULL}. If log dose values are provided the base of
#'   the logarithm should be specified (e.g., \code{exp(1)} for natural logarithm, \code{10}
#'   for base 10).
#' @param control a list of arguments controlling constrained optimisation, maximum iterations,
#'   relative tolerance, and warnings. See \code{\link{drmc}}.
#' @param lowerl a numeric vector of lower limits for all parameters in the model (the default
#'   corresponds to minus infinity for all parameters).
#' @param upperl a numeric vector of upper limits for all parameters in the model (the default
#'   corresponds to plus infinity for all parameters).
#' @param separate logical value indicating whether curves should be fit separately
#'   (independent of each other).
#' @param pshifts a matrix of constants to be added to the matrix of parameters. Default is no
#'   shift for all parameters.
#' @param varcov an optional user-defined known variance-covariance matrix for the responses.
#'   Default is the identity matrix (\code{NULL}), corresponding to independent response values
#'   with a common standard deviation estimated from the data.
#'
#' @return An object of (S3) class \code{"drc"}.
#'
#' @details This function relies on \code{\link{optim}} for minimisation of the negative log
#'   likelihood function. For a continuous response this reduces to least squares estimation.
#'   Response values are assumed to be mutually independent unless \code{varcov} is specified.
#'   For robust estimation MAD (median absolute deviance) is used to estimate the residual
#'   variance. Setting \code{lowerl} and/or \code{upperl} automatically invokes constrained
#'   optimisation. Control arguments may be specified using \code{\link{drmc}}.
#'
#' @seealso \code{\link{drmc}}, \code{\link{LL.4}}, \code{\link{getMeanFunctions}}
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model to the ryegrass data
#' model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' summary(model)
#'
#' @author Christian Ritz and Jens C. Streibig
#'
#' @keywords models nonlinear
"drm" <- function(
formula, curveid, pmodels, weights, data = NULL, subset, fct, 
type = c("continuous", "binomial", "Poisson", "negbin1", "negbin2", "event", "ssd"), bcVal = NULL, bcAdd = 0, 
start, na.action = na.omit, robust = "mean", logDose = NULL, 
control = drmc(), lowerl = NULL, upperl = NULL, separate = FALSE,
pshifts = NULL, varcov = NULL)
{
    ## Matching argument values
    type <- match.arg(type)

    ## Setting na.action option
    op1 <- options(na.action = deparse(substitute(na.action)))
    on.exit(options(op1), add=TRUE)
    
    ## Setting control parameters
    useD <- control$"useD"
    constrained <- control$"constr"
    maxIt <- control$"maxIt"
    optMethod <- control$"method"
    relTol <- control$"relTol"
    warnVal <- control$"warnVal"
    rmNA <- control$"rmNA"
    errorMessage <- control$"errorm"
    noMessage <- control$"noMessage"
    dscaleThres <- control$"dscaleThres"
    rscaleThres <- control$"rscaleThres"
    conCheck <- control$"conCheck"
            
    ## Setting warnings policy
    op2 <- options(warn = warnVal)
    on.exit(options(op2), add=TRUE)

    ## Handling 'start' argument
    if (missing(start)) {selfStart <- TRUE} else {selfStart <- FALSE}

    ## Handling 'fct' argument
    if ( (!is.list(fct)) && (!is.function(fct)) ) {stop("No function or list given in argument 'fct'")}
    if (is.function(fct)) 
    {
        fct <- fct2list(fct, 2)
    }
    
    ## Converting a user specified list
    if (is.null(names(fct))) {fct$"fct" <- fct[[1]]; fct$"ssfct" <- fct[[2]]; fct$"names" <- fct[[3]]}
    
    if (!is.function(fct$"fct")) 
    {
        stop("First entry in list to 'fct' NOT a function")
    } else {
        drcFct <- fct$"fct"
    }
    
    if (is.null(fct$"ssfct")) {noSSfct <- TRUE} else {noSSfct <- FALSE}
    if ((!is.function(fct$"ssfct")) && selfStart)
    {
        stop("Neither self starter function nor starting values provided")
    } else {
        ssfct <- fct$"ssfct"
    }
    
    if (is.null(fct$"names") || (!is.character(fct$"names"))) 
    {
        stop("Parameter names (as vector a strings) are NOT supplied")
    } else {
        parNames <- fct$"names" 
        numNames <- length(parNames)
    }

    ## Checking whether or not first derivates are supplied    
    isDF <- is.function(fct$"deriv1")
    if ( (useD) && (isDF) )
    {
        dfct1 <- fct$"deriv1"
    } else {
        dfct1 <- NULL
    }

    ## Checking whether or not second derivates are supplied    
    if ( (useD) && (is.function(fct$"deriv2")) )
    {
        dfct2 <- fct$"deriv2"
    } else {
        dfct2 <- NULL
    }    

    ## Storing call details
    callDetail <- match.call()

    ## Handling the 'formula', 'curveid' and 'data' arguments
    anName <- deparse(substitute(curveid))  # storing name for later use
    if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)
    if (nchar(anName) < 1) {anName <- "1"}  # in case only one curve is analysed


    mf <- match.call(expand.dots = FALSE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curveid", "data", "subset", "na.action", "weights"), nmf, 0) 

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())  #, globalenv())
    mt <- attr(mf, "terms")

    varNames <- names(mf)[c(2, 1)]  
    varNames0 <- names(mf) 
    # only used once, but mf is overwritten later on
    
    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
    xDim <- ncol(as.matrix(dose))
    resp <- model.response(mf, "numeric")
    if (is.null(resp)) 
    {
        if (xDim > 1) {doseForResp <- dose[, 1]} else {doseForResp <- dose}
        resp <- ppoints(doseForResp, 0.5)[order(doseForResp)]  # just one option
        varNames[1] <- varNames[2] 
        varNames[2] <- "proportion"
    }
    origDose <- dose
    origResp <- resp  # in case of transformation of the response    
    numObs <- length(resp)

    ## Retrieving weights
    wVec <- model.weights(mf)
    if (is.null(wVec))
    {
        wVec <- rep(1, numObs)
    }
    
    ## Finding indices for missing values
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}

    ## Handling "curveid" argument
    assayNo <- model.extract(mf, "curveid") 
    if (is.null(assayNo))  # in case not supplied
    {
        assayNo <- rep(1, numObs)
    }
    uniqueNames <- unique(assayNo)
    colOrder <- order(uniqueNames)
    uniqueNames <- as.character(uniqueNames)

    ## Re-enumerating the levels in 'assayNo' and 'pmodels'
    assayNoOld <- assayNo

    ## Detecting control measurements 
    
    ## Defining helper function     
    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    } 
    assayNo <- colConvert(assayNoOld)
    assayNames <- as.character(unique(assayNoOld))
    numAss <- length(assayNames)
      
#    lenDose <- unlist(lapply(tapply(dose, assayNoOld, unique), length))
    if (xDim > 1) {tempDoseVec <- dose[, 1]} else {tempDoseVec <- dose} 
    uniqueDose <- lapply(tapply(tempDoseVec, assayNoOld, unique), length)
    udNames <- names(uniqueDose[uniqueDose == 1])
    if ( (conCheck) && (length(udNames) > 0) ) 
    {
        cm <- udNames
        if (!noMessage) 
        {
            cat(paste("Control measurements detected for level: ", udNames, "\n", sep = ""))
        }
        
        if (separate)
        {
            stop("Having a common control when fitting separate models does not make sense!\n")
        }
        
        conInd <- assayNoOld %in% udNames
        assayNo[conInd] <- (assayNo[!conInd])[1]
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
        
        ## Updating names, number of curves and the enumeration (starting from 1)
        assayNames <- as.character(unique(assayNoOld[!conInd]))
        numAss <- length(assayNames)
        assayNo <- colConvert(assayNo)
        cm <- NULL
    } else {
        cm <- NULL
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
    }
    
    ## Pooling data from different curves
    if ((separate) && (numAss < 2))
    {
        warning("Only one level: separate = TRUE has no effect", call. = FALSE)
        separate <- FALSE 
    }    
    if ((separate) && (!missing(pmodels)))
    {
        warning("Separate fitting switched off", call. = FALSE)
        separate <- FALSE
    }
    if (separate)
    {
        return(idrm(dose, resp, assayNoOld, wVec, fct, type, control))
    }    
    
    ## Handling "pmodels" argument
    pmodelsList <- list()
    if (missing(pmodels)) 
    {
        if (length(unique(assayNo)) == 1) 
        {
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- matrix(1, numObs, 1)
            }
        } else {
            modelMat <- model.matrix(~ factor(assayNo) - 1, level = unique(assayNo))  # no intercept term
            colnames(modelMat) <- assayNames
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- modelMat
            }
        }   
    } else {
        ## Handling a list or data.frame argument of "pmodels"
        if (is.null(data)) 
        {
            pmodels <- eval(substitute(pmodels), envir = .GlobalEnv)
        } else {
            pmodels <- eval(substitute(pmodels), envir = data, enclos = parent.frame())
        }
   
        if (is.data.frame(pmodels))
        {
            lenCol <- ncol(pmodels)
            pmodelsMat <- matrix(0, numObs, lenCol)    
    
            for (i in 1:lenCol) 
            {
                if (length(unique(pmodels[,i])) == 1) 
                {
                    pmodelsList[[i]] <- matrix(1, numObs, 1)
                    pmodelsMat[,i] <- rep(1, numObs)    
                } else {
                    mf <- eval(model.frame(~factor(pmodels[,i]) - 1), parent.frame())  # converting to factors
                    mt <- attr(mf, "terms")    
    
                    mf2 <- model.matrix(mt, mf)
                    ncmf2 <- ncol(mf2)

                    mf3 <- removeMI(mf2)
                    pmodelsList[[i]] <- mf3
                    pmodelsMat[, i] <- mf3 %*% c(1:ncmf2)
                }
            }
        } else {

            if (is.list(pmodels))
            {   
                lenCol <- length(pmodels)
                pmodelsMat <- matrix(0, length(resp), lenCol)
    
                for (i in 1:lenCol) 
                {
                    if (paste(as.character(pmodels[[i]]), collapse = "") == "~1") 
                    {
                        pmodelsList[[i]] <- matrix(1, numObs, 1)
                        pmodelsMat[,i] <- rep(1, numObs)
                    } else {
                        mf <- eval(model.frame(pmodels[[i]], data=data), parent.frame())   
                        mt <- attr(mf, "terms")    
                        
                        mf2 <- model.matrix(mt, mf)
                        ncmf2 <- ncol(mf2)

                        mf3 <- removeMI(mf2)                    
                        pmodelsList[[i]] <- mf3  
                    
                        pmodelsMat[,i] <- mf3%*%c(1:ncmf2)                    
                    }
                }
            }
        }     
    }

    
    ## Re-setting na.action
    op3 <- options(na.action = "na.omit")  # the default
    on.exit(options(op3), add=TRUE)

    ## Transforming dose value if they are provided as log dose
    if ( !is.null(logDose) && is.numeric(logDose) ) 
    {
       origDose <- dose
       dose <- logDose^dose
    }

    ## Constructing pmodelsList2 from pmodelsList
    pmodelsList2 <- list()
    for (i in 1:numNames)
    {
        if (ncol(pmodelsList[[i]]) > numAss) 
        {
            pmodelsList2[[i]] <- model.matrix(~factor(assayNo) - 1)
            colnames(pmodelsList2[[i]]) <- assayNames 
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])
        }
    }
    
    ## Constructing vectors 'ncclVec' and 'parmPos' used below
    ncclVec <- rep(0, numNames)
    for (i in 1:numNames)
    {
        ncclVec[i] <- ncol(pmodelsList2[[i]])
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])

    ## Constructing parameter names                       
    pnList <- drmParNames(numNames, parNames, pmodelsList2)
    parmVec <- pnList[[1]]
    parmVecA <- pnList[[2]]
    parmVecB <- pnList[[3]]   

    ## Defining with indices for the individual parameters in the model
    parmIndex <- list()
    for (i in 1:numNames)
    {
        parmIndex[[i]] <- parmPos[i] + 1:ncclVec[i]
    }

    ## Scaling of dose and response values 
    scaleFct <- fct$"scaleFct"
    if (!is.null(scaleFct))
    {
        # Defining scaling for dose and response values 
        doseScaling <- 10^(floor(log10(median(dose))))   
        if ( (is.na(doseScaling)) || (doseScaling < dscaleThres) )
        {
            doseScaling <- 1
        }

        respScaling <- 10^(floor(log10(median(resp)))) 
        if ( (is.na(respScaling)) || (respScaling < rscaleThres) || (!identical(type, "continuous")) || (!is.null(bcVal)) )
        {
            respScaling <- 1
        }   

        # Retrieving scaling vector
        longScaleVec <- rep(scaleFct(doseScaling, respScaling), 
                            as.vector(unlist(lapply(parmIndex, length))))
        
    } else {
        doseScaling <- 1
        respScaling <- 1
        
        longScaleVec <- 1
    }      

    ## Constructing vector of initial parameter values
    startVecList <- list()
    
    ## Calculating initial estimates for the parameters using the self starter
    if(!noSSfct)
    {
        startMat <- matrix(0, numAss, numNames)
        lenASS <- length(formals(ssfct))
        if (lenASS > 1)  
        # in case doseScaling and respScaling arguments are available
        # scaling is done inside ssfct()
        {
            doseresp <- data.frame(x = dose, y = origResp)
            ssFct <- function(dframe){ssfct(dframe, doseScaling, respScaling)}
        } else {
        # scaling is explicitly applied to the dose and response values
            doseresp <- data.frame(x = dose / doseScaling, y = origResp / respScaling)
            ssFct <- ssfct
        }

        if (identical(type, "event"))
        {
            dr2 <- doseresp[, 3]
            isFinite <- is.finite(doseresp[, 2])
            respVec <- rep(NA, length(dr2))
            respVec[isFinite] <- cumsum(dr2[isFinite]) / sum(dr2)
            doseresp <- (data.frame(x = doseresp[, 1], y = respVec))[isFinite, ]
        } else {
            isFinite <- is.finite(doseresp[, 2])
        }
        
        ## Finding starting values for each curve
        for (i in 1:numAss)
        {
            indexT1 <- (assayNo[isFinite] == i)
            if (any(indexT1)) 
            {
                logVec <- indexT1
                startMat[i, ] <- ssFct(doseresp[logVec, ])
            } else {
                 startMat[i, ] <- rep(NA, numNames)
            }
    
            ## Identifying a dose response curve only consisting of control measurements
            if (sum(!is.na(startMat[i, ])) == 1) 
            {
                upperPos <- (1:numNames)[!is.na(startMat[i, ])]
            }
        }

        ## Transforming matrix of starting values into a vector
        nrsm <- nrow(startMat)
        for (i in 1:numNames)
        {
            sv <- rep(0, max(nrsm, ncclVec[i]))
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]
            
            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)
            
            startVecList[[i]] <- sv[indVec]
        }
        startVec <- unlist(startVecList)        
    } else {
        startVec <- start
    }
    
    ## Checking the number of start values provided
    if (!selfStart && !noSSfct) 
    {
        lenReq <- length(startVec)  # generated from self starter
        if (length(start) == lenReq) 
        {
            startVec <- start / longScaleVec
        } else {
            stop(paste("Wrong number of initial parameter values. ", lenReq, " values should be supplied", sep = ""))
        }
    }

    ## Converting parameters
    if (selfStart)
    {
        startVec <- drmConvertParm(startVec, startMat, assayNo, pmodelsList2) 
    }

    startVecSc <- startVec

    ## Defining function which converts parameter vector to parameter matrix            
    parmMatrix <- matrix(0, numObs, numNames)
    parm2mat <- function(parm)
    {
        for (i in 1:numNames)
        {
           parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmIndex[[i]]]
        }
        return(parmMatrix)
    }        

    ## Defining model function 
    multCurves <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct", 
                                doseScaling, respScaling, isFinite = rep(TRUE, numObs), pshifts)

    ## Defining first derivative (if available) ... used once in drmEMls()
    if (!is.null(dfct1))
    {
        dmatfct <- function(dose, parm)
        {
            dfct1(dose, parm2mat(parm))
        }
    } else {
        dmatfct <-NULL
    } 

    ## Box-Cox transformation is applied
    if (!is.null(bcVal))
    {
        ## Defining Box-Cox transformation function
        bcfct <- function(x, lambda, bctol, add = bcAdd)
        {
            if (abs(lambda) > bctol)
            {
                return(((x + add)^lambda - 1)/lambda)
            } else {
                return(log(x + add))    
            }
        }
        
        ## Setting the tolerance for Box-Cox transformation being the logarithm transformation 
        ##  (same as in boxcox.default in MASS package)
        bcTol <- 0.02 
        
        resp <- bcfct(resp, bcVal, bcTol)

        multCurves2 <- function(dose, parm)
        {
            bcfct(multCurves(dose, parm), bcVal, bcTol)
        } 
    } else {multCurves2 <- multCurves}

    ## Defining estimation method
    robustFct <- drmRobust(robust, callDetail, numObs, length(startVec))  

    if (type == "continuous")
    {
        ## Ordinary least squares estimation
        estMethod <- drmEMls(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA, dmf = dmatfct, 
        doseScaling = doseScaling, respScaling = respScaling, varcov = varcov)
    }
    if (identical(type, "binomial"))
    {
        estMethod <- drmEMbinomial(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA, 
        doseScaling = doseScaling)        
    } 
    if (identical(type, "Poisson"))
    {
        estMethod <- drmEMPoisson(dose, resp, multCurves2, startVecSc, weightsVec = wVec, 
                                  doseScaling = doseScaling)
    }
    if (identical(type, "negbin1") || identical(type, "negbin2"))
    {
        estMethod <- drmEMnegbin(dose, resp, multCurves2, startVecSc, weightsVec = wVec, 
                                 doseScaling = doseScaling, 
                                 dist.type = ifelse(type == "negbin1", 1, 2))
    }    
    if (identical(type, "event"))
    {
        estMethod <- drmEMeventtime(dose, resp, multCurves2, doseScaling = doseScaling)
    }
    if (identical(type, "ssd"))
    {
        doseScaling <- 1  # dose is the response!
        respScaling <- 1  # no response variable
        longScaleVec <- rep(1, length(longScaleVec))
        multCurves2loc <- modelFunction(dose, parm2mat, fct$"derivx", cm, assayNoOld, upperPos, 
                                        retFct = fct[["retFctDx"]],
                                        doseScaling = doseScaling, respScaling = respScaling, 
                                        isFinite = rep(TRUE, numObs), pshifts)   
        estMethod <- drmEMssd(dose, resp, multCurves2loc, doseScaling = doseScaling, multCurves2 = multCurves2)
    }  
            
    opfct <- estMethod$opfct            

    ## Defining lower and upper limits of parameters
    if (!is.null(lowerl)) 
    {
        if (!is.numeric(lowerl) || !((length(lowerl) == sum(ncclVec)) || (length(lowerl) == numNames)))
        {
            stop("Not correct 'lowerl' argument")
        } else {
            if (length(lowerl) == numNames) 
            {
                lowerLimits <- rep(lowerl, ncclVec)
            } else {
                lowerLimits <- lowerl
            }
        }
        constrained <- TRUE 
        
    } else {  ## In case lower limits are not specified
        lowerLimits <- rep(-Inf, length(startVec))
    }

    if (!is.null(upperl)) 
    {
        if (!is.numeric(upperl) || !((length(upperl) == sum(ncclVec)) || (length(upperl) == numNames)))
        {
            stop("Not correct 'upperl' argument")
        } else {
            if (length(upperl) == numNames) 
            {
                upperLimits <- rep(upperl, ncclVec)
            } else {
                upperLimits <- upperl
            }
        } 
        constrained <- TRUE
                
    } else {  ## In case upper limits are not specified
        upperLimits <- rep(Inf, length(startVec))
    }
    
    lowerLimits <- lowerLimits  / longScaleVec
    upperLimits <- upperLimits  / longScaleVec

    ## Optimisation
    
    ## Setting derivatives
    opdfctTemp <- estMethod$"opdfct1"
    appFct <- function(x, y){tapply(x, y, sum)}   
    
    if (!is.null(opdfctTemp))
    {
        opdfct1 <- function(parm)
        {
            as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo))
        }
    } else {
        opdfct1 <- NULL
    }  

    ## Updating starting values
    startVecSc <- as.vector(startVecSc)  # removing names
    if (identical(type, "negbin1") || identical(type, "negbin2"))
    {
      startVecSc <- c(startVecSc, 1)
      parmVec <- c(parmVec, "O:(Intercept)")
      parmVecA <- c(parmVecA, "O")
      parmVecB <- c(parmVecB, "(Intercept)")
    }  

    ## Optimising the objective function previously defined
    nlsFit <- drmOpt(opfct, opdfct1, startVecSc, optMethod, constrained, warnVal, 
    upperLimits, lowerLimits, errorMessage, maxIt, relTol, parmVec = parmVec, traceVal = control$"trace",
    matchCall = callDetail, silentVal = control$"otrace") 
        
    if (!nlsFit$convergence) {return(nlsFit)}
    
    ## Manipulating after optimisation   
    if (identical(type, "negbin1") || identical(type, "negbin2"))
    {
      longScaleVec <- c(longScaleVec, 1)
    }  
    
    if (identical(type, "event"))
    {
        assayNo0 <- assayNo[isFinite]
        dose0 <- dose[, 2]
        dose1 <- dose0[isFinite]
        dose <- as.vector(unlist(tapply(dose1, assayNo0, function(x){unique(sort(x))})))
        
        ## Rescaling per curve id
        idList <- split(data.frame(dose0, resp), assayNo)
        
        respFct <- function(idListElt)
        {
            doseVec <- idListElt[, 1]
            dose2 <- unique(sort(doseVec))
            orderDose <- order(doseVec)
            resp1 <- tapply(idListElt[orderDose, 2], doseVec[orderDose], sum)  # obtaining one count per time interval
            resp2 <- cumsum(resp1) / sum(resp1)
            
            cbind(dose2, resp2)[is.finite(dose2), , drop = FALSE]
        }
        drList <- lapply(idList, respFct)
        lapList <- lapply(drList, function(x){x[, 1]})
        dose <- as.vector(unlist(lapList))
        resp <- as.vector(unlist(lapply(drList, function(x){x[, 2]})))
        
        splitFactor <- factor(assayNo, exclude = NULL)        
        listCI <- split(splitFactor, splitFactor)
        lenVec <- as.vector(unlist(lapply(lapList, length)))
        plotid <- as.factor(as.vector(unlist(mapply(function(x,y){x[1:y]}, listCI, lenVec))))
        levels(plotid) <- unique(assayNoOld)
    } else {
        plotid <- NULL
    }

    
    if (identical(type, "ssd"))
    {
        if (ncol(as.matrix(dose)) > 1)
        {
            dose2 <- dose[, 2]
            ifDose2 <- is.finite(dose2)
            dose <- dose2[ifDose2]
            resp <- resp[ifDose2]
        } else {
            ifDose2 <- is.finite(dose)  # in case of no censoring  
          
        }
    }
    
    ## Adjusting for pre-fit scaling 
    if (!is.null(scaleFct))
    {
        # Scaling the sums of squares value back
        nlsFit$value <- nlsFit$value * (respScaling^2)
   
        # Scaling estimates and Hessian back
        nlsFit$par <- nlsFit$par * longScaleVec
        nlsFit$hessian <- nlsFit$hessian * (1/outer(longScaleVec/respScaling, longScaleVec/respScaling))
    }
    
    if (!is.null(fct$"retFct"))
    {
        drcFct <- fct$"retFct"(1, 1)  # resetting the scaling
        drcFct1 <- function(dose, parm)
        {
            drcFct(dose, parm2mat(parm)[isFinite, , drop = FALSE])
        }
    }

    # Testing against the ANOVA (F-test)
    nlsSS <- nlsFit$value
    nlsDF <- numObs - length(startVec)

    ## Constructing a plot function
        
    ## Picking parameter estimates for each curve. Does only work for factors not changing within a curve!
    if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames==cm)]} else {iVec <- 1:numAss}
    
    pickCurve <- rep(0, length(iVec))
    for (i in iVec)
    {
       pickCurve[i] <- (1:numObs)[assayNo == i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)

    fixedParm <- (estMethod$"parmfct")(nlsFit)
    parmMat[iVec, ] <- (parm2mat(fixedParm))[pickCurve, ]

    indexMat2 <- parm2mat(1:length(fixedParm))
    indexMat2 <- indexMat2[!duplicated(indexMat2), ]
        
    if (!is.null(cm))
    {
        parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]  
        # 1: simply picking the first row
    }
    rownames(parmMat) <- assayNames


    pmFct <- function(fixedParm)
    {
        if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames == cm)]} else {iVec <- 1:numAss}
    
        if (!is.null(cm))
        {
            parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]  
            # 1: simply picking the first row
        }
        rownames(parmMat) <- assayNames
    
        return(parmMat)
    }
    parmMat <- pmFct(fixedParm)

    ## Defining the plot function    
    pfFct <- function(parmMat)
    {
        plotFct <- function(dose)
        {
            if (is.vector(dose)) 
            {
                lenPts <- length(dose)
            } else {
                lenPts <- nrow(dose)
            }

            curvePts <- matrix(NA, lenPts, ciOrigLength)
            for (i in 1:numAss)
            {
                if (i %in% iVec)
                {
                    parmChosen <- parmMat[i, complete.cases(parmMat[i, ])]  # removing NAs 
                    
                    parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow = TRUE)
                    curvePts[, ciOrigIndex[i]] <- drcFct(dose, parmMat2)
                } else { curvePts[, i] <- rep(NA, lenPts)}
            }
            return(curvePts)
        }
    
        return(plotFct)    
    }
    plotFct <- pfFct(parmMat)


    ## Computation of fitted values and residuals
    if (identical(type, "event") || identical(type, "ssd"))
    {
        multCurves2 <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct", 
                                     doseScaling, respScaling, isFinite)
    }
    predVec <- multCurves2(dose, fixedParm)        
    resVec <- resp - predVec
    resVec[is.nan(predVec)] <- 0    

    diagMat <- matrix(c(predVec, resVec), length(dose), 2)
    colnames(diagMat) <- c("Predicted values", "Residuals")


    ## Adjusting for robust estimation: MAD based on residuals, centered at 0, is used as scale estimate    
    if (robust%in%c("median", "trimmed", "tukey", "winsor"))
    {
        nlsFit$value <- (mad(resVec, 0)^2)*nlsDF 
    }
    if (robust%in%c("lms", "lts"))  # p. 202 i Rousseeuw and Leroy: Robust Regression and Outlier Detection
    {  
        scaleEst <- 1.4826*(1+5/(numObs-length(nlsFit$par)))*sqrt(median(resVec^2))                                 
        w <- (resVec/scaleEst < 2.5)
        nlsFit$value <- sum(w*resVec^2)/(sum(w)-length(nlsFit$par))    
    }
    
    ## Adding meaningful names for robust methods
    robust <- switch(robust, median="median", trimmed="metric trimming", tukey="Tukey's biweight", 
                             winsor="metric Winsorizing", lms="least median of squares",
                             lts="least trimmed squares")


    ## Collecting summary output
    sumVec <- c(NA, NA, NA, nlsSS, nlsDF, numObs)
    sumList <- list(lenData = numObs, 
    alternative = NULL,
    df.residual = numObs - length(startVec))

    ## The data set
    if (!is.null(logDose)) 
    {
        dose <- origDose
    }
    dataSet <- data.frame(origDose, origResp, assayNo, assayNoOld, wVec)
    
    if (identical(type, "event"))
    {
        names(dataSet) <- c(varNames0[c(2, 3, 1)], anName, paste("orig.", anName, sep = ""), "weights")
    } else {
        names(dataSet) <- c(varNames0[c(2, 1)], anName, paste("orig.", anName, sep = ""), "weights")
    }

    ## Matrix of first derivatives evaluated at the parameter estimates
    if (isDF)
    {
        deriv1Mat <- fct$"deriv1"(dose, (parmMat[assayNo, , drop = FALSE])[isFinite, , drop = FALSE])
    } else {
        deriv1Mat <- NULL
    }

    ## Box-Cox information
    if (!is.null(bcVal))
    {
        bcVec <- list(lambda = bcVal, ci = c(NA, NA), bcAdd = bcAdd)
    } else {
        bcVec <- NULL
    }

    ## Parameter estimates
    coefVec <- nlsFit$par
    names(coefVec) <- parmVec
    
    ## Constructing the index matrix
    indexMat <- apply(t(parmMat), 2, function(x){match(x, coefVec)})

    ## Constructing data list
    wName <- callDetail[["weights"]]
    if (is.null(wName)) 
    {
        wName <- "weights"
    } else {
        wName <- deparse(wName)
    }
    dataList <- list(dose = origDose, origResp = as.vector(origResp), weights = wVec, 
    curveid = assayNoOld, resp = as.vector(resp),
    names = list(dName = varNames[1], orName = varNames[2], wName = wName, cNames = anName, rName = ""))
    if (identical(type, "event"))
    {
        dataList <- list(dose = dose, origResp = resp, weights = wVec[isFinite], 
        curveid = assayNoOld[isFinite], plotid = plotid, resp = resp,
        names = list(dName = varNames[1], orName = varNames[2], wName = wName, cNames = anName, rName = ""))
    }
    if (identical(type, "ssd"))
    {
        dataList <- list(dose = dose, origResp = resp, weights = wVec[ifDose2], 
                         curveid = assayNoOld[ifDose2], resp = resp,
                         names = list(dName = varNames[1], orName = varNames[2], 
                                      wName = wName, cNames = anName, rName = ""))
    }

    ## Returning the fit
    returnList <- list(NULL, nlsFit, list(plotFct, logDose), sumVec, startVecSc * longScaleVec, 
    list(parmVec, parmVecA, parmVecB), 
    diagMat, callDetail, dataSet, t(parmMat), fct, robust, estMethod, numObs - length(startVec), 
    sumList, NULL, pmFct, pfFct, type, indexMat, logDose, cm, deriv1Mat, 
    anName, data, wVec, 
    dataList,
    coefVec, bcVec,
    indexMat2)
    
    names(returnList) <- c("varParm", "fit", "curve", "summary", "start", "parNames", 
                           "predres", "call", "data", 
    "parmMat", "fct", "robust", "estMethod", "df.residual", 
    "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat", "logDose", "cm", "deriv1",
    "curveVarNam", "origData", "weights",
    "dataList", "coefficients", "boxcox", "indexMat2")
    class(returnList) <- c("drc")

    return(returnList)
}
