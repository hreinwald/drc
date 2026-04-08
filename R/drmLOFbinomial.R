#' @title Lack-of-fit test for binomial response
#' @keywords internal
"drmLOFbinomial" <- function()
{
    ## Defining a goodness-of-fit test
    gofTest <- function(resp, weights, fitted, dfres)
    {
        ## Removing 0s and 1s in fitted values
        zeroTol <- 1e-12  # no global constant
        indVec <- ( (fitted < zeroTol) | (fitted > 1-zeroTol) )
        dfReduc <- sum(indVec)

        total <- weights
        success <- resp*weights
        expected <- total*fitted

        ## Pearson's statistic (sum of squared Pearson residuals)
        c( sum( ((success - expected)^2 / (expected*(1 - fitted)))[!indVec] ), dfres - dfReduc)
    }


    ## Defining goodness-of-fit function
    anovaTest <- function(formula, ds)
    {
        anovaFit <- glm(formula, family=binomial(link = "logit"), data=ds)
        if (df.residual(anovaFit)>0)
        {
            return(list(test = "lr", anovaFit = anovaFit))
        } else {
            return(NULL)
        }
    }
    anovaTest <- NULL  # lack-of-fit test not meaningful in most situations

    return(list(anovaTest = anovaTest, gofTest = gofTest))
}
