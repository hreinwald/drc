#' Lack-of-fit test for the mean structure based on cumulated residuals
#'
#' The function provides a lack-of-fit test for the mean structure based on cumulated
#' residuals from the model fit.
#'
#' The function provides a graphical model checking of the mean structure in a dose-response
#' model. The graphical display is supplemented by a p-value based on a supremum-type test.
#'
#' The test is applicable even in cases where data are non-normal or exhibit variance
#' heterogeneity.
#'
#' @param object object of class 'drc'.
#' @param noksSim numeric specifying the number of simulations used to obtain the p-value.
#' @param seed numeric specifying the seed value for the random number generator.
#' @param plotit logical indicating whether or not the observed cumulated residual process
#'   should be plotted. Default is to plot the process.
#' @param log character string which should contain \code{"x"} if the x axis is to be
#'   logarithmic, \code{"y"} if the y axis is to be logarithmic and \code{"xy"} or
#'   \code{"yx"} if both axes are to be logarithmic. The empty string \code{""} yields
#'   the original axes.
#' @param bp numeric value specifying the break point below which the dose is zero.
#' @param xlab character string specifying an optional label for the x axis.
#' @param ylab character string specifying an optional label for the y axis.
#' @param ylim numeric vector of length two, containing the lower and upper limit for the y axis.
#' @param ... additional arguments to be passed further to the basic \code{\link{plot}} method.
#'
#' @return A p-value for test of the null hypothesis that the mean structure is appropriate.
#'   Ritz and Martinussen (2009) provide the details.
#'
#' @references Ritz, C and Martinussen, T. (2009) Lack-of-fit tests for assessing mean
#'   structures for continuous dose-response data, \emph{Submitted manuscript}
#'
#' @author Christian Ritz
#'
#' @seealso Other available lack-of-fit tests are the Neill test (\code{\link{neill.test}})
#'   and ANOVA-based test (\code{\link{modelFit}}).
#'
#' @examples
#' ## Fitting a log-logistic model to the dataset 'etmotc'
#' etmotc.m1<-drm(rgr1~dose1, data=etmotc[1:15,], fct=LL.4())
#'
#' ## Test based on cumulated residuals
#' lin.test(etmotc.m1, 1000)
#'
#' @keywords models nonlinear
"lin.test" <- function(object, noksSim = 20, seed = 20070325, plotit = TRUE, log = "", bp = 1e-2, 
xlab, ylab, ylim, ...)
{
    ## Setting seed
    if (!is.null(seed)) {set.seed(seed)}

    ## Retrieving relevant quantities from model fit
    parVec <- coef(object)
    lenpv <- length(parVec)    
    resVec <- residuals(object)
    resVar <- summary(object)$resVar
    xVec <- object$"data"[, 1]
    
    ## Sorting x values
    orderX <- order(xVec)
    xVec <- xVec[orderX]
    resVec <- resVec[orderX]    
    
    lenxv <- length(xVec)
    noObs <- lenxv 
    derVec <- object$"fct"$"deriv1"(xVec, matrix(parVec, lenxv, lenpv, byrow = TRUE))

    ## Calculating W tilde
    eta <- matrix(0, lenxv, lenpv)
    for (i in 1:lenpv)
    {
        eta[, i] <- cumsum(derVec[, i])  # /lenxv
    }
#    term2 <- eta %*% vcov(object) %*% t(matrix(resVec/resVar, 1, lenxv) %*% derVec)

    ## Adjusting in case of replicates
    lenuxv <- length(unique(xVec))
    if (lenuxv < lenxv)
    {
        tempeta <- matrix(0, lenuxv, lenpv)
        for (i in 1:lenpv)
        {
            tempeta[, i] <- as.vector(unlist(tapply(eta[, i], xVec, tail, 1)))
        }
        eta <- tempeta
        lenxv <- lenuxv
        repAdjust <- TRUE
        uxVec <- unique(xVec)
    }  else {
        repAdjust <- FALSE
    }

    ## Generating realisations from the limit distribution
    maxVec <- rep(0, noksSim)
    tempMat <- eta %*% vcov(object)         
    wtMat <- matrix(0, lenxv, noksSim)
    for (i in 1:noksSim)
    {
        rnVec <- rnorm(noObs)
#        wtMat[, i] <- (cumsum(resVec * rnVec) - term2 * rnVec)/sqrt(noObs)
        if (repAdjust)
        {
            term1 <- as.vector(unlist(tapply(cumsum(resVec * rnVec), xVec, tail, 1)))
        } else {
            term1 <- cumsum(resVec * rnVec)
        }
#        print(dim(tempMat))
#        print(dim(matrix(resVec*rnVec/resVar, 1, noObs)))
#        print(dim(derVec))
        term2 <- tempMat %*% t(matrix(resVec*rnVec/resVar, 1, noObs) %*% derVec)
        wti <- (term1 - term2)/sqrt(noObs) 
        wtMat[, i] <- wti
        maxVec[i] <- max(abs(wti))
    }

    ## Calculating the observed process
    if (repAdjust)
    {
        cumRes <- as.vector(unlist(tapply(cumsum(resVec), xVec, tail, 1)))/sqrt(noObs)
        xVec <- unique(xVec)
    } else {
        cumRes <- cumsum(resVec)/sqrt(noObs)
    }

    ## Creating graphical display
    if (plotit)
    {
        ## Finding ranges for x and y in plot
        maxY <- max(c(max(cumRes), max(apply(wtMat, 1, max))))
        minY <- min(c(min(cumRes), min(apply(wtMat, 1, min))))
        
        ## Setting range right in case of log scale
        if (identical(log, "x"))
        {
            lowLim <- bp
            llInd <- xVec > bp
            yStart <- sum(!llInd)
        } else {
            lowLim <- 0
            llInd <- rep(TRUE, length(xVec))
            yStart <- 1
        }
        
        ## Plotting observed and simulated processes
        if (missing(ylim))
        {
            yLim <- c(minY, maxY) 
        } else {
            yLim <- ylim
        }
        
        plot(c(lowLim, xVec[llInd]), c(cumRes[yStart], cumRes[llInd]), type = "s", 
        ylim = yLim, lwd = 2,
        xlab = ifelse(missing(xlab), names(object$"data")[1], xlab), 
        ylab = ifelse(missing(ylab), "Cumulative residuals", ylab), log = log, ...)
        for (i in 1:noksSim)
        {
            lines(c(lowLim, xVec[llInd]), c(wtMat[yStart, i], wtMat[llInd, i]), lty = 2, type = "s")
        }
    }
    
    ## Returning p value for KS test
    sum(max(abs(cumRes)) < maxVec)/noksSim
}

