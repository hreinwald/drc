#' @title The four-parameter Weibull type 1 model
#'
#' @description
#' The general Weibull type 1 model for fitting dose-response data.
#'
#' @details
#' The four-parameter Weibull type 1 model is given by the expression
#' \deqn{f(x) = c + (d - c) \exp(-\exp(b(\log(x) - \log(e))))}
#'
#' The model is sometimes also called the Gompertz model.
#'
#' @param fixed numeric vector of length 4. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 4 giving the names of the
#'   parameters \code{b}, \code{c}, \code{d}, and \code{e}.
#' @param method character string indicating the self starter function to use
#'   (\code{"1"}, \code{"2"}, \code{"3"}, or \code{"4"}).
#' @param ssfct a self starter function to be used. If \code{NULL} (default),
#'   the built-in self starter is used.
#' @param fctName optional character string used internally for the function
#'   name.
#' @param fctText optional character string used internally for the function
#'   text description.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @author Christian Ritz
#'
#' @references
#' Seber, G. A. F. and Wild, C. J. (1989)
#' \emph{Nonlinear Regression}, New York: Wiley & Sons (pp. 338--339).
#'
#' @seealso \code{\link{W1.2}}, \code{\link{W1.3}}, \code{\link{W1.4}},
#'   \code{\link{weibull2}}
#'
#' @keywords models nonlinear
"weibull1" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    
    
#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
    
        parmMat[, 2] + (parmMat[, 3] - parmMat[, 2]) * exp( -exp(parmMat[, 1] *(log(dose) - log(parmMat[, 4]))))
    }


#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos = conPos, fct = confct2))
#    }


#    ## Defining flag to indicate if more general ANOVA model
##    anovaYes <- list(bin = !any(is.na(fixed[c(2,3)])) , cont = TRUE)
#    binVar <- all(fixed[c(2, 3)]==c(0, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}    
#    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining the self starter function
if (FALSE)
{
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[, 1]
        resp3 <- dataFra[, 2]

        startVal <- rep(0, numParm)

        if (is.na(fixed[2]))
        {
            startVal[2] <- min(resp3)  # the lower bound
        } else {
            startVal[2] <- fixed[2]
        }
        
        if (is.na(fixed[3]))
        {
            startVal[3] <- max(resp3)  # the upper bound
        } else {
            startVal[3] <- fixed[3]
        }
        
        
        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA))[notFixed])}  
        # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

#        loglogTrans <- log(-log((resp3-startVal[2] + 0.001)/(startVal[3]-startVal[2])))  # 0.001 to avoid 0 as argument to log
        
#        loglogTrans <- log(-log(abs(resp3 - startVal[2] - startVal[3]/((pi*pi)^2))/(startVal[3] - startVal[2])))
        loglogTrans <- log(-log((resp3 - startVal[2])/(startVal[3] - startVal[2])))

        isFin <- is.finite(loglogTrans)
        loglogTrans <- loglogTrans[isFin]
        dose3 <- dose3[isFin]

#        print(resp3)
#        print(loglogTrans)
#        print(log(dose3))
        loglogFit <- lm(loglogTrans ~ log(dose3))
        
        if (is.na(fixed[4]))
        {
            startVal[4] <- exp(-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter
        } else {
            startVal[4] <- fixed[4]
        }       
#        startVal[4] <- exp(-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter

        if (is.na(fixed[1]))
        {
            startVal[1] <- coef(loglogFit)[2]  # the b parameter
        } else {
            startVal[1] <- fixed[1]
        }       
#        startVal[1] <- coef(loglogFit)[2]  # the b parameter


        ## Avoiding 0 as start value for lower limit (convergence will fail)
        if ( startVal[2] < 1e-12 ) {startVal[2] <- startVal[3]/10}

#        print(startVal)  
        return(startVal[notFixed])
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- weibull1.ssf(method, fixed)
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
    ## Defining derivatives
    deriv1 <- function(dose, parm)
              {
                  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
                  parmMat[, notFixed] <- parm

                  t1 <- parmMat[, 3] - parmMat[, 2]
                  t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
                  t3 <- exp(-t2)

                  derMat <- as.matrix(cbind( -t1 * divAtInf(xlogx(dose/parmMat[, 4], parmMat[, 1]), exp(t2)), 
                                             1 - t3, 
                                             t3, 
                                             t1 * divAtInf(t2, exp(t2)) * parmMat[, 1]/parmMat[, 4] ))
                  return(derMat[, notFixed])
              }
    deriv2 <- NULL

    
    ## Defining the first derivative (in x=dose)
    ##  based on deriv(~c+(d-c)*(exp(-exp(b*(log(x)-log(e))))), "x", function(x, b,c,d,e){})
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm

      .expr1 <- parmMat[, 3] - parmMat[, 2]  # d - c
      .expr6 <- exp(parmMat[, 1] * (log(x) - log(parmMat[, 4])))
      .expr8 <- exp(-.expr6)
      .value <- parmMat[, 2] + .expr1 * .expr8
      .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
      .grad[, "x"] <- -(.expr1 * (.expr8 * (.expr6 * (parmMat[, 1] * (1/x)))))
      .grad
    }
    
#    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)  # function(parm, p, reference, type, ...)
    {        
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)
        
#        if (type == "absolute") {p <- 100*((parmVec[3] - p)/(parmVec[3] - parmVec[2]))}
#        if ( (parmVec[1] < 0) && (reference == "control") ) {p <- 100 - p}
    
        tempVal <- log(-log((100-p)/100))
        EDp <- exp(tempVal/parmVec[1] + log(parmVec[4]))

        EDder <- EDp*c(-tempVal/(parmVec[1]^2), 0, 0, 1/parmVec[4])
    
        return(list(EDp, EDder[notFixed]))
    }

#
#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair)
#    {
#        parmVec1[notFixed] <- parm1
#        parmVec2[notFixed] <- parm2
#
#        tempVal1 <- log(-log((100-pair[1])/100))
#        tempVal2 <- log(-log((100-pair[2])/100))
#    
#        SIpair <- exp(tempVal1/parmVec1[1] + log(parmVec1[4]))/exp(tempVal2/parmVec2[1] + log(parmVec2[4]))
#    
#        SIder1 <- SIpair*c(-tempVal1/(parmVec1[1]*parmVec1[1]), 0, 0, 1/parmVec1[4])
#        SIder2 <- SIpair*c(tempVal2/(parmVec2[1]*parmVec2[1]), 0, 0, -1/parmVec2[4])
#    
#        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
#    }
    

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct, 
#    lowerc=lowerLimits, upperc=upperLimits, confct=confct, anovaYes=anovaYes, 
#    scaleInd = scaleInd,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Weibull (type 1)", fctText),     
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "Weibull-1"
    invisible(returnList)
}


#' @title Two-parameter Weibull type 1 model
#'
#' @description
#' A two-parameter Weibull type 1 model with the lower limit fixed at 0
#' and the upper limit fixed at a specified value (default 1).
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = upper \exp(-\exp(b(\log(x) - \log(e))))}
#'
#' This is mostly used for binomial/quantal responses.
#'
#' @param upper numeric value giving the fixed upper limit. The default is 1.
#' @param fixed numeric vector of length 2. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 2 giving the names of the
#'   parameters. The default is \code{c("b", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{weibull1}}, \code{\link{W1.3}}, \code{\link{W1.4}}
#'
#' @examples
#' earthworms.m1 <- drm(number/total ~ dose, weights = total,
#'   data = earthworms, fct = W1.2(), type = "binomial")
#'
#' @keywords models nonlinear
"W1.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull1(fixed = c(fixed[1], 0, upper, fixed[2]), names = c(names[1], "c", "d", names[2]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Weibull (type 1)", upper), ...))
}

#' @rdname W1.2
w2 <- W1.2

#' @title Three-parameter Weibull type 1 model
#'
#' @description
#' A three-parameter Weibull type 1 model with the lower limit fixed at 0.
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = d \exp(-\exp(b(\log(x) - \log(e))))}
#'
#' This is a special case of the four-parameter Weibull type 1 model
#' where the lower limit is fixed at 0.
#'
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 3 giving the names of the
#'   parameters. The default is \code{c("b", "d", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{weibull1}}, \code{\link{W1.2}}, \code{\link{W1.4}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
#'
#' @keywords models nonlinear
"W1.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull1(fixed = c(fixed[1], 0, fixed[2:3]), names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 1)"), ...))
}

#' @rdname W1.3
w3 <- W1.3

#' @title Three-parameter Weibull type 1 model with upper limit fixed
#'
#' @description
#' A three-parameter Weibull type 1 model with the upper limit fixed
#' (default 1).
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = c + (upper - c) \exp(-\exp(b(\log(x) - \log(e))))}
#'
#' This is a special case of the four-parameter Weibull type 1 model
#' where the upper limit is fixed at a specified value.
#'
#' @param upper numeric value giving the fixed upper limit. The default is 1.
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 3 giving the names of the
#'   parameters. The default is \code{c("b", "c", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @seealso \code{\link{weibull1}}, \code{\link{W1.3}}, \code{\link{W1.4}}
#'
#' @keywords models nonlinear
"W1.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull1(fixed = c(fixed[1:2], upper, fixed[3]), 
    names = c(names[1:2], "d", names[3]), 
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Weibull (type 1)", upper), ...))
}

#' @title Four-parameter Weibull type 1 model
#'
#' @description
#' A four-parameter Weibull type 1 model.
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = c + (d - c) \exp(-\exp(b(\log(x) - \log(e))))}
#'
#' @param fixed numeric vector of length 4. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 4 giving the names of the
#'   parameters. The default is \code{c("b", "c", "d", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @references
#' Seber, G. A. F. and Wild, C. J. (1989)
#' \emph{Nonlinear Regression}, New York: Wiley & Sons (pp. 338--339).
#'
#' Ritz, C. (2009)
#' Towards a unified approach to dose-response modeling in ecotoxicology.
#' \emph{Environ Toxicol Chem}, \bold{29}, 220--229.
#'
#' @seealso \code{\link{weibull1}}, \code{\link{W1.2}}, \code{\link{W1.3}}
#'
#' @examples
#' terbuthylazin.m1 <- drm(rgr ~ dose, data = terbuthylazin, fct = W1.4())
#'
#' @keywords models nonlinear
"W1.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}

    return(weibull1(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]),
    fctText = "Weibull (type 1)", ...))
}

#' @rdname W1.4
w4 <- W1.4


#' @title Two-parameter exponential decay model
#'
#' @description
#' A two-parameter exponential decay model with the slope parameter \code{b}
#' fixed at 1 and the lower limit fixed at 0.
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = d \exp(-x/e)}
#'
#' This is a special case of the Weibull type 1 model
#' (\code{\link{weibull1}}) with the slope fixed at 1 and the lower limit
#' fixed at 0.
#'
#' @param fixed numeric vector of length 2. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 2 giving the names of the
#'   parameters. The default is \code{c("d", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @references
#' Seber, G. A. F. and Wild, C. J. (1989)
#' \emph{Nonlinear Regression}, New York: Wiley & Sons (pp. 338--339).
#'
#' @seealso \code{\link{EXD.3}}, \code{\link{AR.2}}, \code{\link{AR.3}},
#'   \code{\link{weibull1}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = EXD.2())
#'
#' @keywords models nonlinear
"EXD.2" <-
function(fixed = c(NA, NA), names = c("d", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull1(fixed = c(1, 0, fixed[1:2]), 
    names = c("b", "c", names[1:2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Exponential decay"), ...))
}

#' @title Three-parameter exponential decay model
#'
#' @description
#' A three-parameter exponential decay model with the slope parameter \code{b}
#' fixed at 1.
#'
#' @details
#' The model is given by the expression
#' \deqn{f(x) = c + (d - c) \exp(-x/e)}
#'
#' This is a special case of the Weibull type 1 model
#' (\code{\link{weibull1}}) with the slope fixed at 1.
#'
#' @param fixed numeric vector of length 3. Specifies which parameters are
#'   fixed and at what value. Use \code{NA} for parameters that are not fixed.
#' @param names character vector of length 3 giving the names of the
#'   parameters. The default is \code{c("c", "d", "e")}.
#' @param \dots additional arguments passed to \code{\link{weibull1}}.
#'
#' @return A list of class \code{Weibull-1} containing the nonlinear function,
#'   self starter function, and parameter names.
#'
#' @references
#' Seber, G. A. F. and Wild, C. J. (1989)
#' \emph{Nonlinear Regression}, New York: Wiley & Sons (pp. 338--339).
#'
#' @seealso \code{\link{EXD.2}}, \code{\link{AR.2}}, \code{\link{AR.3}},
#'   \code{\link{weibull1}}
#'
#' @examples
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = EXD.3())
#'
#' @keywords models nonlinear
"EXD.3" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}

    return(weibull1(fixed = c(1, fixed[1:3]), 
    names = c("b", names[1:3]),
    fctName = as.character(match.call()[[1]]),
    fctText = "Shifted exponential decay", ...))
}
