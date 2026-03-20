"drmEMbinomial" <- 
function(dose, resp, multCurves, startVec, robustFct, weights, rmNA, zeroTol = 1e-12, 
doseScaling = 1, respScaling = 1)
{
    ## Finding indices for doses that give contribution to likelihood function
    iv <- ( (multCurves(dose/doseScaling, startVec) > zeroTol) & (multCurves(dose/doseScaling, startVec) < 1-zeroTol) )

    ## Defining the objective function                
#' @title EM algorithm for binomial response
#' @keywords internal
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
        prob <- multCurves(dose / doseScaling, c)
        omZT <- 1 - zeroTol
        prob[prob > omZT] <- omZT
        prob[prob < zeroTol] <- zeroTol
        -sum((resp * weights) * log(prob / (1 - prob)) + (weights * log(1 - prob)))
    }    
    
    ## Defining self starter function
    ssfct <- NULL

    ## Defining the log likelihood function
    llfct <- function(object)
    {
        total <- (object$"data")[, 5]
        success <- total*(object$"data")[, 2]    
        
        c(sum(log(choose(total, success))) - object$"fit"$"ovalue",
        object$"sumList"$"lenData" - df.residual(object))
    }
       
    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }

    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, 
    vcovfct = vcovfct, parmfct = parmfct))
}


