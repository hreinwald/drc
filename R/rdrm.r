"rdrm" <- function(nosim, fct, mpar, xerror, xpar = 1, yerror = "rnorm", ypar = c(0, 1), 
onlyY = FALSE)
{        
    ## Constructing the predictor values
    if (is.numeric(xerror))
    {
        x <- xerror
    } else {
        xFun <- match.fun(xerror)
        x <- do.call(xFun, as.list(xpar))
    }
    lenx <- length(x)
    x <- sort(x)
    x <- rep(x, nosim)
    xMat <- matrix(x, nosim, lenx, byrow = TRUE)
    
    ## Constructing the mean dose-response
    meanVec <- fct$fct(x, matrix(mpar, lenx*nosim, length(mpar), byrow = TRUE))
    
    ## Constructing the simulated response values
    yFun <- match.fun(yerror)
    if (yerror == "rbinom")
    {
        if (length(ypar) == 1)
        {
            ypar <- rep(ypar, lenx*nosim)
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        } else {
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        }
        errorVec <- yFun(lenx*nosim, ypar, meanVec)
        
        yMat <- matrix(errorVec, nosim, lenx, byrow = TRUE)

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {    
            return(list(x = xMat, w = wMat, y = yMat))
        }
    }  else {
        errorVec <- do.call(yFun, c(list(lenx*nosim), as.list(ypar)))
        
        yMat <- matrix(meanVec, nosim, lenx, byrow = TRUE) + errorVec

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {
            return(list(x = xMat, y = yMat))
        }
    }
}
