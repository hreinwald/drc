#' Creating isobolograms
#'
#' \code{isobole} displays isobole based on EC/ED50 estimates from a log-logistic model.
#' Additionally isoboles determined by the concentration addition model, Hewlett's model
#' and Voelund's model can be added to the plot.
#'
#' The model fits to be supplied as first and optionally second argument are obtained
#' using \code{\link{mixture}} and \code{\link{drm}}.
#'
#' @param object1 object of class 'drc' where EC/ED50 parameters vary freely.
#' @param object2 object of class 'drc' where EC/ED50 parameters vary according to Hewlett's model.
#' @param exchange numeric. The exchange rate between the two substances.
#' @param cifactor numeric. The factor to be used in the confidence intervals. Default is 2,
#'   but 1 has been used in publications.
#' @param ename character string. The name of the EC/ED50 variable.
#' @param xaxis character string. Is the mixture "0:100" or "100:0" on the x axis?
#' @param xlab an optional label for the x axis.
#' @param ylab an optional label for the y axis.
#' @param xlim a numeric vector of length two, containing the lower and upper limit for the x axis.
#' @param ylim a numeric vector of length two, containing the lower and upper limit for the y axis.
#' @param ... Additional graphical parameters.
#'
#' @return No value is returned. Only used for the side effect: the isobologram shown.
#'
#' @references Ritz, C. and Streibig, J. C. (2014) From additivity to synergism - A
#'   modelling perspective \emph{Synergy}, \bold{1}, 22--29.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"isobole" <- 
function(object1, object2, exchange = 1, cifactor = 2, ename = "e", xaxis = "100",  
xlab, ylab, xlim, ylim, ...)
{

    parmVec <- coef(object1)
    namesPV <- names(parmVec)
#    lenNPV <- length(namesPV)

    indVec <- regexpr(paste(ename, ":", sep = ""), namesPV, fixed = TRUE) > 0
    eVec <- parmVec[indVec]
    seVec <- (summary(object1)$"coefficients")[indVec, 2]

#    edMat <- ED(object1, 50, display = FALSE, multcomp = TRUE)[["EDdisplay"]]
    edMat <- ED(object1, 50, display = FALSE)
    eVec <- (as.vector(edMat[, 1]))  # stripping off names
    seVec <- (as.vector(edMat[, 2]))  # stripping off names

    mixProp <- unique(object1$data[, 4])
    mixProp <- mixProp[ (mixProp >= 0) & (mixProp <= 100) ] / 100  # removing control level
#    print(mixProp)

#    posOnRay <- function(len, slope)
#    {
#        xPos <- sqrt( len^2 / (1+slope^2) )
#        yPos <- slope*xPos
#        
#        yPos[!is.finite(slope)] <- len[!is.finite(slope)]
#        
#        list(xPos, yPos)
#    }

#    if (identical(xaxis, "0")) 
#    {
#        eVec <- rev(eVec)
#        seVec <- rev(seVec)
#        mixProp <- rev(mixProp)
#    }

    Ex <- eVec * mixProp
    Ey <- eVec * (1-mixProp) * exchange
#    print(Ex)
#    print(Ey)

    lowerE <- eVec - cifactor * seVec
#    lowerEx <- lowerE*mixProp
#    lowerEx <- lowerE*cos(mixProp*pi/2)
#    lowerEy <- lowerE*(1-mixProp)*exchange    
#    lowerEy <- lowerE*sin(mixProp*pi/2)*exchange    

    upperE <- eVec + cifactor * seVec
#    upperEx <- upperE*mixProp
#    upperEx <- upperE*cos(mixProp*pi/2)
#    upperEy <- upperE*(1-mixProp)*exchange
#    upperEy <- upperE*sin(mixProp*pi/2)*exchange

#    lowerE <- eVec - 2 * seVec
    lowerEx <- lowerE * mixProp
    lowerEy <- lowerE * (1 - mixProp) * exchange
#    upperE <- eVec + 2 * seVec
    upperEx <- upperE * mixProp
    upperEy <- upperE * (1 - mixProp) * exchange


    ## Defining plotting frame
    if (missing(xlim)) {xLimits <- c(0, max(upperEx))} else {xLimits <- xlim}       
    if (missing(ylim)) {yLimits <- c(0, max(upperEy))} else {yLimits <- ylim}
    
    if (missing(xlab)) {xLab <- ifelse(identical(xaxis, "100"), "100", "0")} else {xLab <- xlab}
    if (missing(ylab)) {yLab <- ifelse(identical(xaxis, "100"), "0", "100")} else {yLab <- ylab}

    ## Swapping axes
    if (identical(xaxis, "0")) 
    {
        ETemp <- Ex
        Ex <- Ey
        Ey <- ETemp
        
        lowerTemp <- lowerEx
        lowerEx <- lowerEy
        lowerEy <- lowerTemp

        upperTemp <- upperEx
        upperEx <- upperEy
        upperEy <- upperTemp

        limTemp <- xLimits
        xLimits <- yLimits
        yLimits <- limTemp 
        
        mixProp <- 1 - mixProp
    }
    
    ## Drawing the plot frame    
    plot(0, type = "n", xlab = xLab, ylab = yLab, xlim = xLimits, ylim = yLimits, ...)  # empty plot
    
    ## Plotting rays in first quadrant
    raySlopes <- (1 - mixProp) / mixProp
#    raySlopes <- mixProp/(1 - mixProp) 
#    raySlopes <- tan(mixProp * pi/2)
    for (i in raySlopes[is.finite(raySlopes)]) {abline(0, exchange*i, lty = 3)}
    abline(v = 0, lty = 3)  # adding vertical line (for infinite slope)       

#    raySlopes <- mixProp/(1 - mixProp)
#    for (i in raySlopes[is.finite(raySlopes)]) 
#    {
#        abline(0, exchange * i, lty = 3)
#    }

    ## Plotting ED50 values with confidence intervals  
    points(Ex, Ey, pch = 19)
    segments(lowerEx, lowerEy, upperEx, upperEy, lwd = 2)    
    
#    points(eVec*cos(mixProp*pi/2), eVec*sin(mixProp*pi/2)*exchange, pch = 19)

#    katx <- function(eVal, slope) {cos(atan(slope))*eVal}
#    katy <- function(eVal, slope) {sin(atan(slope))*eVal}    

#    points(katx(eVec, raySlopes), katy(eVec, raySlopes)*exchange, pch = 19)
#    segments(lowerEx, lowerEy, upperEx, upperEy, lwd = 2)
# old    segments(katx(lowerE, raySlopes), katy(lowerE, raySlopes)*exchange, 
#        katx(upperE, raySlopes), katy(upperE, raySlopes)*exchange, lwd = 2)

    if (!missing(object2))
    {

        ## Retrieving parameter estimates from fit of Hewlett's model
        parmVec <- coef(object2)
        namesPV <- names(parmVec)
#        lenNPV <- length(namesPV)

        curveStr1 <- paste("I(1/(", object1$curveVarNam, "/100))", sep = "")
        curveStr2 <- paste("I(1/(1 - ", object1$curveVarNam, "/100))", sep = "")
        ED50.1 <- parmVec[regexpr(curveStr1, namesPV, fixed = TRUE) > 0]
        ED50.2 <- parmVec[regexpr(curveStr2, namesPV, fixed = TRUE) > 0]
        ## maybe use an extractor like 'EDmix' instead of the above 4 lines
  
        ## Plotting the isobole based on the Voelund model
        if (identical(object2$fct$name, "voelund"))         
        {
            Voelund0 <- function(t0, m0, eta, prop)
            {
                if (prop == 1) 
                {
                    retVal <- log(t0) 
                } else if (prop == 0) 
                {
                    retVal <- log(m0)
                } else {
                    propr <- (1 - prop) / prop
                    tempVal <- (1+t0/m0*propr)^(1-eta[1]) + (t0/m0*propr)^eta[2] * (1+t0/m0*propr)^(1-eta[2])
                    retVal <- log((t0 / tempVal) * (1 + propr))
                } 
                retVal
            }
            Voelund <- Vectorize(Voelund0, "prop")
                                     
            mixVec <- seq(1, 0, length = 100)  # oops hardcoded "100"
            voeVec <- exp(Voelund(ED50.1, ED50.2, tail(parmVec, 2), mixVec))        
            tvals <- mixVec * voeVec
            mvals <- (1 - mixVec) * voeVec

            xVal <- tvals
            yVal <- exchange * mvals

        } else {
            Hewlett <- function(t, t0, m0, lambda)
            {
                retVal <- (m0^(1/lambda) - (t*m0/t0)^(1/lambda))^lambda
                retVal[t > t0] <- 0
                retVal
            }    
            if (identical(object2$fct$"name", "hewlett")) 
            {
                lambda <- tail(parmVec , 1)            
            } 
            if (identical(object2$fct$"name", "ca"))
            {
                lambda <- 1
            }
            xVal <- seq(0, ED50.1, length = 100)
            yVal <- exchange * Hewlett(xVal, ED50.1, ED50.2, lambda)

        }
        if (identical(xaxis, "0")) 
        {
            tempVal <- xVal
            xVal <- yVal
            yVal <- tempVal
        }    
        lines(xVal, yVal, ...)        
    }
#    invisible()
}
