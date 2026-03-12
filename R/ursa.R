#' Universal Response Surface Approach (URSA) for Drug Interaction
#'
#' URSA provides a parametric approach for modelling the joint action of several
#' agents. The model allows quantification of synergistic effects through a single
#' parameter. The model function is defined implicitly through an appropriate equation.
#'
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value
#'   they are fixed. NAs for parameters that are not fixed.
#' @param names a vector of character strings giving the names of the parameters.
#'   The default is reasonable.
#' @param ssfct a self starter function to be used (optional).
#'
#' @return A list containing the nonlinear function, the self starter function,
#'   and the parameter names.
#'
#' @references
#'   Greco, W. R. and Park H. S. and Rustum, Y. M. (1990) Application of a New
#'   Approach for the Quantitation of Drug Synergism to the Combination of
#'   cis-Diamminedichloroplatinum and 1-beta-D-Arabinofuranosylcytosine,
#'   \emph{Cancer Research}, \bold{50}, 5318--5327.
#'
#'   Greco, W. R. Bravo, G. and Parsons, J. C. (1995) The Search for Synergy:
#'   A Critical Review from a Response Surface Perspective,
#'   \emph{Pharmacological Reviews}, \bold{47}, Issue 2, 331--385.
#'
#' @author Christian Ritz after an idea by Hugo Ceulemans.
#'
#' @seealso Other models for fitting mixture data: \code{\link{mixture}}.
#'
#' @examples
#' \dontrun{
#' d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2,
#' 2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 20, 20, 20, 20,
#' 20, 50, 50, 50, 50, 50)
#' d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2,
#' 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2,
#' 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)
#' effect <- c(106, 99.2, 115, 79.2, 70.1, 49, 21, 3.83, 74.2,
#' 71.5, 48.1, 30.9, 16.3, 76.3, 48.8, 44.5, 15.5, 3.21, 56.7,
#' 47.5, 26.8, 16.9, 3.25, 46.7, 35.6, 21.5, 11.1, 2.94, 24.8,
#' 21.6, 17.3, 7.78, 1.84, 13.6, 11.1, 6.43, 3.34, 0.89)
#' greco <- data.frame(d1, d2, effect)
#' greco.m1 <- drm(effect ~ d1 + d2, data = greco,
#'   fct = ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA)))
#' summary(greco.m1)
#' }
#'
#' @keywords models nonlinear
"ursa" <- function(
fixed = rep(NA, 7), names = c("b1", "b2", "c", "d", "e1", "e2", "f"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 7
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    bisec <- function(fu, fuLow, fuHigh)
    {
	      for(k in 1:25) 
        {
	          fuMiddle <- 0.5 * (fuLow + fuHigh)
		        if (fu(fuMiddle) > 0) 
            {
			          fuHigh <- fuMiddle
      		  } else {
			          fuLow <- fuMiddle
      		  }
    	  }    
        list(root = fuMiddle)
    }
    
    ## Defining the nonlinear model function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        parmMat[, 5] <- parmMat[, 5] / dose[, 1]
        parmMat[, 6] <- parmMat[, 6] / dose[, 2]

        applyImplicitFct <- function(parmVec)
        {
#            print(parmVec)
            
            if ((!is.finite(parmVec[5])) && (!is.finite(parmVec[6])))
            {
                return(parmVec[4])
            } else {

                implicitFct <- function(e)
                {
                    scaledEffect <- (e - parmVec[3]) / (parmVec[4] - e)
                    recSlope1 <- 1 / parmVec[1]
                    recSlope2 <- 1 / parmVec[2] 
        
                    1/(parmVec[5] * (scaledEffect^recSlope1)) + 1/(parmVec[6] * (scaledEffect^recSlope2)) + 
                    parmVec[7]/(parmVec[5] * parmVec[6] * (scaledEffect^(recSlope1/2 + recSlope2/2))) - 1
                } 
        
#                print(c(implicitFct(parmVec[3]*1.01), implicitFct(parmVec[4]*0.99), parmVec[4]*0.99))
        
#                reducFactor0 <- max(c((1/parmVec[5])^parmVec[1] + 1, (1/parmVec[6])^parmVec[2] + 1))
#                reducFactor <- max(c(0.99, reducFactor0 / (1 + reducFactor0)))
#                print(c((1/parmVec[5])^parmVec[1], (1/parmVec[6])^parmVec[2], reducFactor0, reducFactor))
#                bisection <- try(uniroot(implicitFct, c(parmVec[3] * (2 - reducFactor), parmVec[4] * reducFactor)), silent = FALSE)                        
#                bisection <- try(uniroot(implicitFct, c(parmVec[3]*1.01, parmVec[4]*0.99)), silent = FALSE)
                bisection <- try(bisec(implicitFct, parmVec[3], parmVec[4]), silent = TRUE)
                if (inherits(bisection, "try-error")) 
                {
                    return(NA)
                } else {
                    return(bisection$"root")
                } 
            }       
        }      

#        print(apply(parmMat, 1, applyImplicitFct))                  
        apply(parmMat, 1, applyImplicitFct)
    }    

    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
#            initval <- c((llogistic()$ssfct(dframe))[c(1, 1:4, 4)], 0.5) * c(-1, -1, rep(1, 5))
            
            startLL.d1 <- as.vector(coef(drm(dframe[, c(3,1)], fct = LL.4())))
            startLL.d2 <- as.vector(coef(drm(dframe[, c(3,2)], fct = LL.4())))

            if (startLL.d1[1] < 0)  # condition in terms of standard "drc" parameter "b"
            {
                initVal <- c(startLL.d1, startLL.d2, 0)[c(1, 5, 3, 2, 4, 8, 9)]          
            } else {
                initVal <- c(startLL.d1, startLL.d2, 0)[c(1, 5, 2, 3, 4, 8, 9)] * c(-1, -1, rep(1, 5))           
            }
    
            return(initVal[notFixed])              
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
#    scaleFct <- function(doseScaling, respScaling)
#    {        
#        c(1, 1, respScaling, respScaling, doseScaling, doseScaling, 1)[notFixed]
#    }

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct = edfct, sifct = sifct,  # scaleFct = scaleFct,
    name = "ursa",
    text = "URSA", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "ursa"
    invisible(returnList)
}
