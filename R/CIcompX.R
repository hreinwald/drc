## Calculating combination indices for x and y axes

#' Calculation of combination index for binary mixtures
#'
#' For single mixture data, combination indices for effective doses as well as effects
#' may be calculated. This is an extended version of \code{\link{CIcomp}}.
#'
#' @param mixProp a numeric value between 0 and 1 specifying the mixture proportion/ratio.
#' @param modelList a list containing 3 model fits using \code{\link{drm}}: the mixture model fit
#'   first, followed by the 2 pure substance model fits.
#' @param EDvec a numeric vector of effect levels (percentages between 0 and 100).
#' @param EDonly logical. If TRUE, only combination indices for effective doses are calculated.
#'
#' @return A list with components \code{Effx}, \code{Effy} (unless \code{EDonly = TRUE}),
#'   \code{CAx}, \code{CAy} (unless \code{EDonly = TRUE}), and \code{EDvec}.
#'
#' @references Martin-Betancor, K. and Ritz, C. and Fernandez-Pinas, F. and Leganes, F. and
#'   Rodea-Palomares, I. (2015) Defining an additivity framework for mixture research in
#'   inducible whole-cell biosensors, \emph{Scientific Reports} \bold{17200}.
#'
#' @author Christian Ritz and Ismael Rodea-Palomares
#'
#' @seealso \code{\link{CIcomp}}, \code{\link{plotFACI}}, \code{\link{mixture}}
#'
#' @keywords models nonlinear
#' @concept antagonism mixture synergy
CIcompX <- function(mixProp, modelList, EDvec, EDonly = FALSE)
{
    ## Checking the input
    if ( (mixProp < 0) | (mixProp > 1) ) {stop("Mixture proportion should be between 0 and 1")}
    if ( (!is.list(modelList)) | (length(modelList) != 3) ) {stop("Exactly 3 model fits should be provided in a list")}
    if (length(EDvec) < 1) {stop("At least effective dose level should be specified")} 

    ## Estimating effective doses and effects
    ese12Vec <- ED(modelList[[1]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]  # the mixture
    ese1Vec <- ED(modelList[[2]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]
    ese2Vec <- ED(modelList[[3]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]    
    eseMat <- as.matrix(cbind(ese12Vec[, 1], ese1Vec[, 1], ese2Vec[, 1], ese12Vec[, 2], ese1Vec[, 2], ese2Vec[, 2]))
    rownames(eseMat) <- as.character(EDvec)
    colnames(eseMat) <- c("ED.mix", "ED1", "ED2", "SE.mix", "SE1", "SE2")

    pred12 <- predict(modelList[[1]], data.frame(ese12Vec[, 1]), se.fit = TRUE)  
    # Note: Ignoring uncertainty in the ED values!
    pred1 <- predict(modelList[[2]], data.frame(ese1Vec[, 1]), se.fit = TRUE)
    pred2 <- predict(modelList[[3]], data.frame(ese2Vec[, 1]), se.fit = TRUE) 
    
    ## In case only a single ED level is specified 
    ## (as predict() then returns a vector, not a matrix)
    if (!is.matrix(pred12)) {pred12 <- matrix(pred12, nrow = 1)}
    if (!is.matrix(pred1)) {pred1 <- matrix(pred1, nrow = 1)}
    if (!is.matrix(pred2)) {pred2 <- matrix(pred2, nrow = 1)}
       
    predMat <- as.matrix(cbind(pred12[, 1], pred1[, 1], pred2[, 1], 
                               pred12[, 2], pred1[, 2], pred2[, 2]))
    rownames(predMat) <- as.character(EDvec)
    colnames(predMat) <- c("E.mix", "E1", "E2", "SE.mix", "SE1", "SE2")

    ## Calculating combination index for effective doses
    xCAfct <- function(eseVec)
    {
        combInd <- (mixProp*eseVec[1]/eseVec[2]) + ((1-mixProp)*eseVec[1]/eseVec[3]) 
        caDiff <- combInd - 1

        derivFct <- function(ecVec)
        {
            derivComp <- c(
            mixProp / ecVec[2] + (1-mixProp) / ecVec[3], 
            -mixProp * ecVec[1] / (ecVec[2]^2), 
            -(1 - mixProp) * ecVec[1] / ((ecVec[3]^2)))
            
            derivComp
        }
        derivVec <- derivFct(eseVec[1:3])
        diagVec <- diag(eseVec[4:6]^2)
 #       seCI <- sqrt(as.vector((-derivVec[2:3]) %*% (diagVec[2:3, 2:3]) %*% (-derivVec[2:3])))
        seDiff <- sqrt(as.vector(derivVec %*% diagVec %*% derivVec))

        derivFct2 <- function(ecVec)
        {
            addEff <- (mixProp / ecVec[2] + (1-mixProp) / ecVec[3])^{-1}
            derivComp2 <- c(0, (addEff^2) * mixProp / (ecVec[2]^2), (addEff^2) * (1 - mixProp) / (ecVec[3]^2))
        
            c(addEff, derivComp2)
        }
        dfRes2 <- derivFct2(eseVec[1:3])
        derivVec2 <- dfRes2[2:4]
        seDiff2 <- sqrt(as.vector(derivVec2 %*% diagVec %*% derivVec2))
#        print(seDiff2)

        retVec <- c(combInd, seDiff, c(combInd - 1.96 * seDiff, combInd + 1.96 * seDiff), 
                    caDiff, 2 * (1 - pnorm(abs(caDiff / seDiff))), dfRes2[1], seDiff2)
        names(retVec) <- c("combInd", "SE",  
                           "lowCI", "highCI", 
                           "CAdiff", "CAdiffp", "PredAdd", "sePredAdd")
        retVec
    }

# CI = {1/ [(pi/EX1) + (1-pi)/Ex2]}/eseVec[1]
# CI = {1/ [(pi/eseVec[2]) + (1-pi)/eseVec[3]]}/eseVec[1]
# CI = {1/ [(pi*eseVec[1]/eseVec[2]) + (1-pi)*eseVec[1]/eseVec[3]]}    
    
    ## Calculating combination index for effects
    yCAfct <- function(eseVec)
    {
        denomi <- (mixProp*eseVec[1]/eseVec[2]) + ((1-mixProp)*eseVec[1]/eseVec[3])
        combInd <- 1 / denomi
        caDiff <- combInd - 1

        derivFct <- function(ecVec)
        {
            derivComp <- (-1/(denomi^2)) * c(
            mixProp / ecVec[2] + (1-mixProp) / ecVec[3] , 
            -mixProp * ecVec[1] / (ecVec[2]^2), 
            -(1 - mixProp) * ecVec[1] / (ecVec[3]^2))
        }
        derivVec <- derivFct(eseVec[1:3])
        diagVec <- diag(eseVec[4:6]^2)
 #       seCI <- sqrt(as.vector((-derivVec[2:3]) %*% (diagVec[2:3, 2:3]) %*% (-derivVec[2:3])))
        seDiff <- sqrt(as.vector(derivVec %*% diagVec %*% derivVec))

        derivFct2 <- function(ecVec)
        {
            addEff <- (mixProp / ecVec[2] + (1-mixProp) / ecVec[3])^{-1}
            derivComp2 <- c(0, (addEff^2) * mixProp / (ecVec[2]^2), (addEff^2) * (1 - mixProp) / (ecVec[3]^2))
        
            c(addEff, derivComp2)
        }
        dfRes2 <- derivFct2(eseVec[1:3])
        derivVec2 <- dfRes2[2:4]
        seDiff2 <- sqrt(as.vector(derivVec2 %*% diagVec %*% derivVec2))                
        
        retVec <- c(combInd, seDiff, 
                    c(combInd - 1.96 * seDiff, combInd + 1.96 * seDiff),
                    caDiff, 2 * (1 - pnorm(abs(caDiff / seDiff))), dfRes2[1], seDiff2)
        names(retVec) <- c("combInd", "SE",
                           "lowCI", "highCI", 
                           "CAdiff", "CAdiffp", "PredAdd", "sePredAdd")
        retVec
    }
    CAxMat <- t(apply(eseMat, 1, xCAfct))
    rownames(CAxMat) <- EDvec    
    CAyMat <- t(apply(predMat, 1, yCAfct))  # not yCAfct
    rownames(CAyMat) <- EDvec   
    
    if (EDonly)
    {
        list(Effx = eseMat, CAx = CAxMat, EDvec = EDvec)
    } else {
        list(Effx = eseMat, Effy = predMat, CAx = CAxMat, CAy = CAyMat, EDvec = EDvec)
    }
}

#' Classical combination index for effective doses
#'
#' Calculates the classical combination index for effective doses in binary mixture experiments.
#'
#' @param mixProp a numeric value between 0 and 1 specifying the mixture proportion/ratio.
#' @param modelList a list containing 3 model fits using \code{\link{drm}}: the mixture model fit
#'   first, followed by the 2 pure substance model fits.
#' @param EDvec a numeric vector of effect levels (percentages between 0 and 100).
#'
#' @return A matrix with one row per ED value. Columns contain estimated combination indices,
#'   their standard errors and 95% confidence intervals, p-value for testing CI=1, estimated
#'   ED values for the mixture data and assuming concentration addition (CA) with corresponding
#'   standard errors.
#'
#' @references Martin-Betancor, K. and Ritz, C. and Fernandez-Pinas, F. and Leganes, F. and
#'   Rodea-Palomares, I. (2015) Defining an additivity framework for mixture research in
#'   inducible whole-cell biosensors, \emph{Scientific Reports} \bold{17200}.
#'
#' @author Christian Ritz and Ismael Rodea-Palomares
#'
#' @seealso \code{\link{CIcompX}}, \code{\link{plotFACI}}, \code{\link{mixture}}
#'
#' @examples
#' ## Fitting marginal models for the 2 pure substances
#' acidiq.0 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 0), fct = LL.4())
#' acidiq.100 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 100), fct = LL.4())
#'
#' ## Fitting model for single mixture with ratio 17:83
#' acidiq.17 <- drm(rgr ~ dose, data = subset(acidiq, pct == 17 | pct == 0), fct = LL.4())
#'
#' ## Calculation of combination indices based on ED10, ED20, ED50
#' CIcomp(0.17, list(acidiq.17, acidiq.0, acidiq.100), c(10, 20, 50))
#'
#' @keywords models nonlinear
#' @concept antagonism mixture synergy
CIcomp <- function(mixProp, modelList, EDvec)
{
    resLst <- CIcompX(mixProp, modelList, EDvec, EDonly = FALSE)
    resMt <- matrix(NA, 5, 5)

    resMt[c(1,2,4), 1:2] <- matrix(resLst[["Effx"]][c(2, 3, 1, 5, 6, 4)], 3, 2)
    resMt[3, 1:2] <- resLst[["CAx"]][7:8]
    resMt[5, ] <- resLst[["CAx"]][c(1:2, 4:6)]
    
    resMt <- cbind(resLst[["CAx"]][, -5], resLst[["Effx"]][, c(1, 4)])
    colnames(resMt)[6:7] <- c("ED.CA", "SE.CA")

    # colnames(resMt) <- c("Est", "SE", "CIlow", "CIupp", "p-val")
    # rownames(resMt) <- c("ED.A", "ED.B", "ED.CApred", "ED.mix", "CombInd")    
    
    resMt
}


#' Plot combination index as a function of fraction affected
#'
#' Visualizes the combination index from \code{\link{CIcompX}} as a function of the fraction affected.
#'
#' @param effList a list as returned by \code{\link{CIcompX}}.
#' @param indAxis character string. Either "ED" for effective doses or "EF" for effects.
#' @param caRef logical. If TRUE (default), a reference line for concentration addition is drawn.
#' @param showPoints logical. If TRUE, estimated combination indices are plotted as points.
#' @param add logical. If TRUE, the plot is added to an existing plot.
#' @param ylim numeric vector of length 2 giving the range for the y axis.
#' @param ... additional graphical arguments.
#'
#' @return Invisibly returns the plot matrix of combination index values.
#'
#' @author Christian Ritz and Ismael Rodea-Palomares
#'
#' @seealso \code{\link{CIcompX}}, \code{\link{CIcomp}}
#'
#' @keywords models nonlinear
plotFACI <- function(effList, indAxis = c("ED", "EF"), caRef = TRUE, 
                     showPoints = FALSE, add = FALSE, ylim, ...)
{
    indAxis <- match.arg(indAxis)
#    indMat <- CIcompX(mixProp, modelList, faValues)
  
    faValues <- effList[["EDvec"]]
    minfa <- min(faValues)
    faValues[faValues < 0] <- -(100 - abs(faValues[faValues < 0])) 
    
#    if (indAxis == "x") {plotMat <- indMat[[1]]} else {plotMat <- indMat[[2]]}
    plotMat <- switch(indAxis, ED = effList[["CAx"]], EF = effList[["CAy"]])
    xVec <- as.numeric(rownames(plotMat))
    xVec[faValues < 0] <- rev(xVec[faValues < 0])
    xVec[faValues > 0] <- rev(xVec[faValues > 0])
    yVec <- plotMat[, 1]
    seVec <- plotMat[, 2]
    if (caRef) 
    {
        yLimits <- c(0, max(yVec))
    } else {
        yLimits <- range(yVec)
    }
    if (!missing(ylim))
    {
        yLimits <- ylim
    }
    
    if (!add)
    {
        plot(xVec, yVec, type = "l", xlim = c(minfa - 1, max(faValues) + 1), ylim = yLimits, 
        xlab = "Fraction affected", ylab = "Combination index", ...)
        
        abline(h = 1, lty = 2, lwd = 2)
    
    } else {
        lines(xVec, yVec, ...)   
    }
    dispersion(xVec, yVec, ulim = 1.96 * seVec, arrow.gap = 0.15, ...)
    if (showPoints) 
    {
        points(xVec, yVec, pch = 1)
    }
    
    invisible(plotMat)
}
