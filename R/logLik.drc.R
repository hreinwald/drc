#' @title Extracting the log likelihood
#'
#' @description
#' \code{logLik} extracts the value of the log likelihood function evaluated
#' at the parameter estimates.
#'
#' @param object an object of class 'drc'.
#' @param ... additional arguments.
#'
#' @return The evaluated log likelihood as a numeric value and the
#'   corresponding degrees of freedom as well as the number of observations
#'   as attributes.
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())
#' logLik(ryegrass.m1)
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"logLik.drc" <- function(object, ...)
{
    ## Retrieving the value of the log likelihood function evaluated at the parameter estimates

    if (inherits(object, "bindrc"))
    {
        llVal <- object$loglik[3]
        attr(llVal, "df") <- object$loglik[5] - object$loglik[4]
        attr(llVal, "nobs") <- nrow(object$data) # <--  extension provided by Tobias Verbeke
    } else {

#        llVal <- summary(object)[[4]][1]
        
#        if (object$"type"=="continuous")
#        {
            loglik <- (object$estMethod$llfct)(object)
            
            llVal <- loglik[1]
            attr(llVal, "df") <- loglik[2]   
            attr(llVal, "nobs") <- nrow(object$data) # <--  extension provided by Tobias Verbeke       
#        }

#        if (object$"type"=="binomial")
#        {
#            degfre <- object$summary[6]
#        
#            total <- (object$"data")[,5]
#            success <- total*(object$"data")[,2]    
#            llVal <- sum(log(choose(total, success))) - object$fit$ofvalue

#            attr(llVal, "df") <- object[[4]][7] - object[[4]][6]            
#        }
        
        
#        numVarPar <- 1  # + 1 to add variance parameter
#        if (object$"type"=="binomial") {numVarPar <- 0}
#        attr(llVal, "df") <- object[[4]][7] - object[[4]][6] + numVarPar
#        attr(llVal, "df") <- object[[4]][6]
    }

    class(llVal) <- "logLik"
    return(llVal)
}

