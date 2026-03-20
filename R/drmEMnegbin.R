#' @importFrom stats dnbinom
"drmEMnegbin" <- 
function(dose, resp, multCurves, startVec, weightsVec, doseScaling = 1, dist.type = 1)
{

    ## Defining the objective function  
    if (dist.type == 1)
    {  
        opfct <- function(cVal)
        {
            sizeVal <- tail(cVal, 1)
            pVal <- 1 / (1 + weightsVec * multCurves(dose / doseScaling, head(cVal, -1)) * exp(sizeVal)) 
            -sum(dnbinom(resp, exp(-sizeVal), pVal, log = TRUE))
        }
    }

    if (dist.type == 2)
    {  
        opfct <- function(cVal)
        {
            sizeVal <- tail(cVal, 1)
            pVal <- 1 / (1 + weightsVec * exp(sizeVal)) 
            -sum(dnbinom(resp, exp(-sizeVal) * multCurves(dose / doseScaling, head(cVal, -1)), 
                         pVal, log = TRUE))
        }
    }
  
    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        c(-object$"fit"$value, object$"sumList"$"df.residual"
        )  # adding scale constant
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
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


#' @title EM algorithm for negative binomial
#' @keywords internal
"drmLOFnegbin" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
