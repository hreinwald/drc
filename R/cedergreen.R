#' @title The Cedergreen-Ritz-Streibig model
#'
#' @description
#' \code{cedergreen} provides a very general way of specifying the Cedergreen-Ritz-Streibig
#' modified log-logistic model for describing hormesis, under various constraints on the parameters.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed.
#'   NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters (should not contain ":").
#'   The order of the parameters is: b, c, d, e, f.
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used.
#' @param alpha numeric value between 0 and 1, reflecting the steepness of the hormesis peak.
#'   This argument must be specified.
#' @param fctName optional character string used internally by convenience functions.
#' @param fctText optional character string used internally by convenience functions.
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = c + \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
#' which is a five-parameter model (alpha is fixed).
#'
#' @return A list containing the non-linear function, the self starter function
#'   and the parameter names.
#'
#' @references
#' Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005)
#' Improved empirical models describing hormesis,
#' \emph{Environmental Toxicology and Chemistry} \bold{24}, 3166--3172.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{CRS.4a}}, \code{\link{CRS.5a}}, \code{\link{ucedergreen}}, \code{\link{drm}}
#'
#' @examples
#' lettuce.crsm1 <- drm(weight ~ conc, data = lettuce, fct = CRS.6())
#' summary(lettuce.crsm1)
#'
#' @keywords models nonlinear
"cedergreen" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL, 
alpha, fctName, fctText)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}

    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2] + (parmMat[,3] - parmMat[,2] + parmMat[,5]*exp(-1/(dose^alpha)))/(1 + exp(parmMat[,1]*(log(dose) - log(parmMat[,4]))))
    }


#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos=conPos, fct=confct2))
#    }
#
#
#    ## Defining flag to indicate if more general ANOVA model
##    anovaYes <- TRUE
#    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}
#    anovaYes <- list(bin = binVar, cont = TRUE)

    ## Defining self starter function
if (FALSE)
{    
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

#        startVal[3]<-max(resp3)+0.001  # the d parameter
#        startVal[2]<-min(resp3)-0.001  # the c parameter
        startVal[3] <- 1.05 * resp3[which.min(dose2)]
        startVal[2] <- 0.95 * min(resp3)
        
#        startVal[!notFixed] <- fixed[!notFixed] 

        if (length(unique(dose2))==1) {return((c(NA,NA,startVal[3],NA,NA))[notFixed])}

        indexT2<-(dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}
        dose3<-dose2[indexT2]
        resp3<-resp3[indexT2]

        logitTrans<-log((startVal[3]-resp3)/(resp3-startVal[2] + 0.001))  # 0.001 to avoid 0 in the denominator
        logitFit<-lm(logitTrans~log(dose3))
        startVal[4]<-exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1]<-coef(logitFit)[2]  # the b parameter

#        startVal[5] <- 0  # the f parameter
        ## Solving equation at x=e
        startVal[5] <- (2*(median(resp3) - startVal[2]) - (startVal[3] - startVal[2]))*exp(1/(startVal[4]^alpha))

        return(startVal[notFixed])
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
#        ssfct <- cedergreen.ssf(method, fixed, alpha)       
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)   
            initval[5] <- (2*(median(dframe[, 2])-initval[2])-(initval[3]-initval[2]))*exp(1/(initval[4]^alpha))
    
            return(initval[notFixed])
        }        
    }
       
    ## Defining names
    names <- names[notFixed]

#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }


    ## Defining derivatives
    
#    ## Constructing a helper function
#    xlogx <- function(x, p)
#    {
#        lv <- (x < 1e-12)
#        nlv <- !lv
#        
#        rv <- rep(0, length(x))
#        
#        xlv <- x[lv] 
#        rv[lv] <- log(xlv^(xlv^p[lv]))
#        
#        xnlv <- x[nlv]
#        rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
#    
#        rv
#    }
    
    ## Specifying the derivatives    
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t0 <- exp(-1/(dose^alpha))
        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*t0
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1 - 1/t3, 
               1/t3, 
               t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               t0/t3 )[, notFixed]
    }
        
    deriv2 <- NULL


    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function    
#    edfct <- function(parm, p, lower = 1e-4, upper = 10000, ...)  # upper2=1000)
    edfct <- function(parm, respl, reference, type, lower = 1e-4, upper = 10000, ...)  # upper2=1000)
    {
#        if (is.null(upper)) {upper <- 1000}
#        if (missing(upper2)) {upper2 <- 1000}
        interval <- c(lower, upper) 
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type, TRUE)  # FALSE)  Changed 2010-06-02 after e-mail from Claire       
        tempVal <- (100-p)/100

        helpFct <- function(dose) {parmVec[2]+(parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))/(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))}
#        doseVec <- exp(seq(-upper2, upper2, length=1000))
        doseVec <- exp(seq(log(interval[1]), log(interval[2]), length=1000))
        maxAt <- doseVec[which.max(helpFct(doseVec))]
#        print(maxAt)
#        print(upper)
    
        eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*exp(-1/(dose^alpha))/(parmVec[3]-parmVec[2]))}
        EDp <- uniroot(eqn, lower=maxAt, upper=upper)$root

        EDdose <- EDp
        tempVal1 <- exp(parmVec[1]*(log(EDdose)-log(parmVec[4])))
        tempVal2 <- parmVec[3]-parmVec[2]
        derParm <- c(tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), -parmVec[5]*exp(-1/(EDdose^alpha))/((tempVal2)^2),
                     parmVec[5]*exp(-1/(EDdose^alpha))/((tempVal2)^2), -tempVal*tempVal1*parmVec[1]/parmVec[4],
                     -exp(-1/(EDdose^alpha))/tempVal2)
        derDose <- tempVal*tempVal1*parmVec[1]/EDdose-parmVec[5]/tempVal2*exp(-1/(EDdose^alpha))/(EDdose^(1+alpha))*alpha 

        EDder <- derParm/derDose
        
        return(list(EDp, EDder[notFixed]))
    }

#
#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair, upper = 10000, interval = c(1e-4, 10000))
#    {
##        if (is.null(upper)) {upper <- 1000}
##        if (missing(upper2)) {upper2 <- 1000}
#    
#        parmVec1[notFixed] <- parm1
#        parmVec2[notFixed] <- parm2
#    
#        tempVal1 <- (100-pair[1])/100
#        tempVal2 <- (100-pair[2])/100
#
##        doseVec <- exp(seq(-upper2, upper2, length=max(c(1000, upper2))))
#        doseVec <- exp(seq(log(interval[1]), log(interval[2]), length=1000))        
#        helpFct1 <- function(dose) 
#                    {
#                        parmVec1[2]+(parmVec1[3]-parmVec1[2]+parmVec1[5]*exp(-1/(dose^alpha)))/(1+exp(parmVec1[1]*(log(dose)-log(parmVec1[4]))))
#                    }
#        maxAt1 <- doseVec[which.max(helpFct1(doseVec))]
#
#        helpFct2 <- function(dose) 
#                    {
#                        parmVec2[2]+(parmVec2[3]-parmVec2[2]+parmVec2[5]*exp(-1/(dose^alpha)))/(1+exp(parmVec2[1]*(log(dose)-log(parmVec2[4]))))
#                    }
#        maxAt2 <- doseVec[which.max(helpFct2(doseVec))]
#
#        eqn1 <- function(dose) {tempVal1*(1+exp(parmVec1[1]*(log(dose)-log(parmVec1[4]))))-(1+parmVec1[5]*exp(-1/(dose^alpha))/(parmVec1[3]-parmVec1[2]))}
#        EDp1 <- uniroot(eqn1, lower=maxAt1, upper=upper)$root
#        eqn2 <- function(dose) {tempVal2*(1+exp(parmVec2[1]*(log(dose)-log(parmVec2[4]))))-(1+parmVec2[5]*exp(-1/(dose^alpha))/(parmVec2[3]-parmVec2[2]))}
#        EDp2 <- uniroot(eqn2, lower=maxAt2, upper=upper)$root
#
#        SIpair <- EDp1/EDp2
#
#        EDdose1 <- EDp1
#        EDdose2 <- EDp2
#        tempVal11 <- exp(parmVec1[1]*(log(EDdose1)-log(parmVec1[4])))
#        tempVal12 <- parmVec1[3]-parmVec1[2]
#        derParm1 <- c(tempVal1*tempVal11*(log(EDdose1)-log(parmVec1[4])), -parmVec1[5]*exp(-1/(EDdose1^alpha))/((tempVal12)^2),
#                     parmVec1[5]*exp(-1/(EDdose1^alpha))/((tempVal12)^2), -tempVal1*tempVal11*parmVec1[1]/parmVec1[4],
#                     -exp(-1/(EDdose1^alpha))/tempVal12)
#        derDose1 <- tempVal1*tempVal11*parmVec1[1]/EDdose1-parmVec1[5]/tempVal12*exp(-1/(EDdose1^alpha))/(EDdose1^(1+alpha))*alpha 
#
#        SIder1 <- (derParm1/derDose1)/EDp2
#
#        tempVal21 <- exp(parmVec2[1]*(log(EDdose2)-log(parmVec2[4])))
#        tempVal22 <- parmVec2[3]-parmVec2[2]
#        derParm2 <- c(tempVal2*tempVal21*(log(EDdose2)-log(parmVec2[4])), -parmVec2[5]*exp(-1/(EDdose2^alpha))/((tempVal22)^2),
#                     parmVec2[5]*exp(-1/(EDdose2^alpha))/((tempVal22)^2), -tempVal2*tempVal21*parmVec2[1]/parmVec2[4],
#                     -exp(-1/(EDdose2^alpha))/tempVal22)
#        derDose2 <- tempVal2*tempVal21*parmVec2[1]/EDdose2-parmVec2[5]/tempVal22*exp(-1/(EDdose2^alpha))/(EDdose2^(1+alpha))*alpha 
#
#        SIder2 <- (derParm2/derDose2)*(-EDp1/(EDp2*EDp2))
#
#        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
#    }


    ## Finding the maximal hormesis
    maxfct <- function(parm, lower = 1e-3, upper = 1000)
    {
#        if (is.null(upper)) {upper <- 1000}
#        if (is.null(interval)) {interval <- c(1e-3, 1000)}            
#        alpha <- 0.5
        parmVec[notFixed] <- parm
        
        optfct <- function(t)
        {
            expTerm1 <- parmVec[5]*exp(-1/(t^alpha))
            expTerm2 <- exp(parmVec[1]*(log(t)-log(parmVec[4])))
            
            return(expTerm1*alpha/(t^(alpha+1))*(1+expTerm2)-(parmVec[3]-parmVec[2]+expTerm1)*expTerm2*parmVec[1]/t)
        }
    
        ED1 <- edfct(parm, 1, lower, upper)[[1]]
               
        doseVec <- exp(seq(log(1e-6), log(ED1), length = 100))
#        print((doseVec[optfct(doseVec)>0])[1])

        maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec)>0])[1], ED1))$root
        return(c(maxDose, fct(maxDose, matrix(parm, 1, length(names)))))
    }

    
    returnList <- 
    list(fct=fct, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2,  # lowerc=lowerLimits, upperc=upperLimits, 
    edfct=edfct, maxfct=maxfct, 
#    scaleInd=scaleInd, anovaYes=anovaYes, confct=confct,
#    name = "cedergreen",
#    text = "Cedergreen-Ritz-Streibig", 
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Cedergreen-Ritz-Streibig", fctText),     
    noParm = sum(is.na(fixed)))

    class(returnList) <- "mllogistic"
    invisible(returnList)
}

#' @title Cedergreen-Ritz-Streibig model with lower limit 0 (alpha=1)
#'
#' @description
#' Four-parameter CRS hormesis model with the lower limit fixed at 0 and alpha=1.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.5a}}, \code{\link{UCRS.4a}}
#'
#' @examples
#' lettuce.crsm1 <- drm(lettuce[,c(2,1)], fct = CRS.4a())
#' summary(lettuce.crsm1)
#' ED(lettuce.crsm1, c(50))
#'
#' @keywords models nonlinear
"CRS.4a" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 1, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1)",
    ...))
}

#' @title Alias for CRS.4a
#' @description \code{ml3a} is an alias for \code{\link{CRS.4a}}.
#' @seealso \code{\link{CRS.4a}}
#' @keywords models nonlinear
ml3a <- CRS.4a

#' @title Cedergreen-Ritz-Streibig model with lower limit 0 (alpha=0.5)
#'
#' @description
#' Four-parameter CRS hormesis model with the lower limit fixed at 0 and alpha=0.5.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.4a}}, \code{\link{CRS.5b}}
#'
#' @keywords models nonlinear
"CRS.4b" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 0.5, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=.5)",
    ...))
}

#' @title Alias for CRS.4b
#' @description \code{ml3b} is an alias for \code{\link{CRS.4b}}.
#' @seealso \code{\link{CRS.4b}}
#' @keywords models nonlinear
ml3b <- CRS.4b

#' @title Cedergreen-Ritz-Streibig model with lower limit 0 (alpha=0.25)
#'
#' @description
#' Four-parameter CRS hormesis model with the lower limit fixed at 0 and alpha=0.25.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.4a}}, \code{\link{CRS.5c}}
#'
#' @keywords models nonlinear
"CRS.4c" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==4)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, 0, NA, NA, NA), names = c(names[1], "c", names[2:4]), alpha = 0.25, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=.25)",
    ...))
}

#' @title Alias for CRS.4c
#' @description \code{ml3c} is an alias for \code{\link{CRS.4c}}.
#' @seealso \code{\link{CRS.4c}}
#' @keywords models nonlinear
ml3c <- CRS.4c

#' @title Cedergreen-Ritz-Streibig five-parameter model (alpha=1)
#'
#' @description
#' Five-parameter CRS hormesis model with alpha=1.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.4a}}, \code{\link{UCRS.5a}}
#'
#' @examples
#' lettuce.m1 <- drm(lettuce[,c(2,1)], fct = CRS.5a())
#' summary(lettuce.m1)
#' ED(lettuce.m1, c(50))
#'
#' @keywords models nonlinear
"CRS.5a" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==5)) {stop("Not correct 'names' argument")}
 
    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 1, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=1)",
    ...))
}

#' @title Alias for CRS.5a
#' @description \code{ml4a} is an alias for \code{\link{CRS.5a}}.
#' @seealso \code{\link{CRS.5a}}
#' @keywords models nonlinear
ml4a <- CRS.5a

#' @title Cedergreen-Ritz-Streibig five-parameter model (alpha=0.5)
#'
#' @description
#' Five-parameter CRS hormesis model with alpha=0.5.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.4b}}, \code{\link{CRS.5a}}
#'
#' @keywords models nonlinear
"CRS.5b" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==5)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 0.5, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=.5)",
    ...))
}

#' @title Alias for CRS.5b
#' @description \code{ml4b} is an alias for \code{\link{CRS.5b}}.
#' @seealso \code{\link{CRS.5b}}
#' @keywords models nonlinear
ml4b <- CRS.5b

#' @title Cedergreen-Ritz-Streibig five-parameter model (alpha=0.25)
#'
#' @description
#' Five-parameter CRS hormesis model with alpha=0.25.
#'
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{cedergreen}}.
#'
#' @return A list (see \code{\link{cedergreen}}).
#'
#' @seealso \code{\link{cedergreen}}, \code{\link{CRS.4c}}, \code{\link{CRS.5a}}
#'
#' @keywords models nonlinear
"CRS.5c" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names)==5)) {stop("Not correct 'names' argument")}

    return(cedergreen(fixed = c(NA, NA, NA, NA, NA), names = names, alpha = 0.25, 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Cedergreen-Ritz-Streibig (alpha=.25)",
    ...))
}

#' @title Alias for CRS.5c
#' @description \code{ml4c} is an alias for \code{\link{CRS.5c}}.
#' @seealso \code{\link{CRS.5c}}
#' @keywords models nonlinear
ml4c <- CRS.5c

