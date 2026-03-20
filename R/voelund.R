#' Voelund Mixture Model
#'
#' Provides the Voelund model for describing the joint action of two compounds
#' in binary mixture experiments. Used internally by \code{\link{mixture}}.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters
#'   (should not contain ":").
#' @param method character string indicating the self starter function to use.
#' @param ssfct a self starter function to be used (optional).
#' @param eps numeric tolerance for handling zero dose values.
#'
#' @return A list containing the nonlinear model function, the self starter function,
#'   and the parameter names.
#'
#' @author Christian Ritz
#'
#' @seealso \code{\link{mixture}}, \code{\link{hewlett}}
#'
#' @keywords internal
"voelund" <- function(
fixed = c(NA, NA, NA, NA, NA, NA, NA), 
names = c("b", "c", "d", "e", "f", "g", "h"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
eps = 1e-10)
{
    ## Checking arguments
    numParm <- 7
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

        ratio<- parmMat[, 4]/parmMat[, 5]
        tmp <- (1+ratio)^(1-parmMat[, 6])+((ratio)^parmMat[, 7])*((1+ratio)^(1-parmMat[, 7]))
        loge <- log(parmMat[, 4]/tmp)        
        
        loge[!is.finite(parmMat[, 4])] <- log(parmMat[!is.finite(parmMat[, 4]), 5])
        loge[!is.finite(parmMat[, 5])] <- log(parmMat[!is.finite(parmMat[, 5]), 4])
        
        retVec <- parmMat[, 2]+(parmMat[, 3]-parmMat[, 2])/(1+exp(parmMat[, 1]*(log(dose)-loge)))
        ## Handling the case dose=0 where "loge" may become NaN due to the mixture encoding (pct in glymet)
        retVec[dose < eps] <- parmMat[dose < eps, 3]
        retVec        
    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- c((llogistic()$ssfct(dframe))[c(1:4, 4)], 1, 1)
    
            return(initval[notFixed])
        }        
    } 
 
   
    ## Defining names
    names <- names[notFixed]

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL

    ## Scale function
    scaleFct <- function(doseScaling, respScaling)
    {        
        c(1, respScaling, respScaling, doseScaling, doseScaling, 1, 1)[notFixed]
    }

    returnList <- 
    list(fct=fct, ssfct=ssfct, names=names, deriv1=deriv1, deriv2=deriv2, 
    edfct=edfct, sifct=sifct, scaleFct=scaleFct,
    name = "voelund",
    text = "Voelund mixture", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "Voelund"
    invisible(returnList)
}
