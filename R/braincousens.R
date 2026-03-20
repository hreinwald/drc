#' @title The Brain-Cousens hormesis models
#'
#' @description
#' \code{braincousens} provides a very general way of specifying Brain-Cousens'
#' modified log-logistic model for describing hormesis, under various constraints on the parameters.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed.
#'   NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters (should not contain ":").
#'   The order of the parameters is: b, c, d, e, f.
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used.
#' @param fctName optional character string used internally by convenience functions.
#' @param fctText optional character string used internally by convenience functions.
#'
#' @details
#' The Brain-Cousens model is given by the expression
#' \deqn{f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
#' which is a five-parameter model.
#'
#' @return A list containing the non-linear function, the self starter function,
#'   the parameter names and additional model specific objects.
#'
#' @references
#' Brain, P. and Cousens, R. (1989) An equation to describe dose responses
#' where there is stimulation of growth at low doses,
#' \emph{Weed Research}, \bold{29}, 93--96.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{BC.4}}, \code{\link{BC.5}}, \code{\link{drm}}
#'
#' @keywords models nonlinear
"braincousens" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    


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
        
        parmMat[,2]+(parmMat[,3]+parmMat[,5]*dose-parmMat[,2])/(1+exp(parmMat[,1]*(log(dose)-log(parmMat[,4]))))
    }


    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)   
            initval[5] <- 0
    
            return(initval[notFixed])
        }           
    }
   
    ## Defining names
    names <- names[notFixed]


    ## Defining derivatives
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*dose
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1 - 1/t3, 
               1/t3, 
               t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               dose/t3 )[, notFixed]
    }
        
    deriv2 <- NULL


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, lower = 1e-3, upper = 1000, ...)
    {
        interval <- c(lower, upper)     
     
        parmVec[notFixed] <- parm

        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- (100-p)/100

        helpEqn <- function(dose) 
        {
            expVal <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            parmVec[5]*(1+expVal*(1-parmVec[1]))-(parmVec[3]-parmVec[2])*expVal*parmVec[1]/dose
        }
        maxAt <- uniroot(helpEqn, interval)$root
    
        eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*dose/(parmVec[3]-parmVec[2]))}
        EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root

        EDdose <- EDp
        tempVal1 <- exp(parmVec[1]*(log(EDdose)-log(parmVec[4])))
        tempVal2 <- parmVec[3]-parmVec[2]
        derParm <- c(tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), -parmVec[5]*EDdose/((tempVal2)^2),
                     parmVec[5]*EDdose/((tempVal2)^2), -tempVal*tempVal1*parmVec[1]/parmVec[4],
                     -EDdose/tempVal2)
        derDose <- tempVal*tempVal1*parmVec[1]/EDdose-parmVec[5]/tempVal2 

        EDder <- derParm/derDose
        
        return(list(EDp, EDder[notFixed]))
    }


    ## Finding the maximal hormesis
    maxfct <- function(parm, lower = 1e-3, upper = 1000)
    {
        parmVec[notFixed] <- parm
        if (parmVec[1]<1) {stop("Brain-Cousens model with b<1 not meaningful")}
        if (parmVec[5]<0) {stop("Brain-Cousens model with f<0 not meaningful")}
        
        optfct <- function(t)
        {
            expTerm1 <- parmVec[5]*t
            expTerm2 <- exp(parmVec[1]*(log(t)-log(parmVec[4])))
            
            return(parmVec[5]*(1+expTerm2)-(parmVec[3]-parmVec[2]+expTerm1)*expTerm2*parmVec[1]/t)
        }
    
        ED1 <- edfct(parm, 1, lower, upper)[[1]]
               
        doseVec <- exp(seq(log(1e-6), log(ED1), length = 100))

        maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec)>0])[1], ED1))$root
        return(c(maxDose, fct(maxDose, matrix(parm, 1, length(names)))))
    }


    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct = edfct, maxfct = maxfct, 
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Brain-Cousens (hormesis)", fctText),
    noParm = sum(is.na(fixed)))

    class(returnList) <- "braincousens"
    invisible(returnList)
}


#' @title Four-parameter Brain-Cousens hormesis model
#'
#' @description
#' \code{BC.4} provides the Brain-Cousens modified log-logistic model with the lower limit fixed at 0.
#'
#' @param fixed numeric vector of length 4 specifying fixed parameters (NAs for free parameters).
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{braincousens}}.
#'
#' @return A list (see \code{\link{braincousens}}).
#'
#' @references
#' van Ewijk, P. H. and Hoekstra, J. A. (1993)
#' Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
#' \emph{Ecotoxicology and Environmental Safety}, \bold{25}, 25--32.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{braincousens}}, \code{\link{BC.5}}
#'
#' @examples
#' lettuce.bcm2 <- drm(weight ~ conc, data = lettuce, fct = BC.4())
#' summary(lettuce.bcm2)
#' ED(lettuce.bcm2, c(50))
#'
#' @keywords models nonlinear
"BC.4" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}
 
    return(braincousens(names=c(names[1], "c", names[2:4]), fixed = c(fixed[1], 0, fixed[2:4]),
    fctName = as.character(match.call()[[1]]),
    fctText = "Brain-Cousens (hormesis) with lower limit fixed at 0", ...))
}

#' @title Alias for BC.4
#' @description \code{bcl3} is an alias for \code{\link{BC.4}}.
#' @param fixed numeric vector of length 4 specifying fixed parameters (NAs for free parameters).
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{braincousens}}.
#' @seealso \code{\link{BC.4}}
#' @keywords models nonlinear
bcl3 <- BC.4

#' @title Five-parameter Brain-Cousens hormesis model
#'
#' @description
#' \code{BC.5} provides the full five-parameter Brain-Cousens modified log-logistic model
#' for describing hormesis.
#'
#' @param fixed numeric vector of length 5 specifying fixed parameters (NAs for free parameters).
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{braincousens}}.
#'
#' @return A list (see \code{\link{braincousens}}).
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{braincousens}}, \code{\link{BC.4}}
#'
#' @examples
#' lettuce.bcm1 <- drm(weight ~ conc, data = lettuce, fct = BC.5())
#' modelFit(lettuce.bcm1)
#' plot(lettuce.bcm1)
#'
#' @keywords models nonlinear
"BC.5" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(braincousens(names = names, fixed = fixed,
    fctName = as.character(match.call()[[1]]), ...))
}

#' @title Alias for BC.5
#' @description \code{bcl4} is an alias for \code{\link{BC.5}}.
#' @param fixed numeric vector of length 5 specifying fixed parameters (NAs for free parameters).
#' @param names a vector of character strings giving the names of the parameters.
#' @param ... additional arguments passed to \code{\link{braincousens}}.
#' @seealso \code{\link{BC.5}}
#' @keywords models nonlinear
bcl4 <- BC.5
