"drmEMPoisson" <- 
function(dose, resp, multCurves, startVec, weightsVec, doseScaling = 1)
{

    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
        lambda <- weightsVec * multCurves(dose / doseScaling, c)
        return( -sum(-lambda + resp*log(lambda)))
    }    

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        c(
        -object$"fit"$value + sum(log(gamma(resp+1))),
        object$"sumList"$"df.residual"
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


#' @title EM algorithm for Poisson response
#' @keywords internal
"drmLOFPoisson" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
